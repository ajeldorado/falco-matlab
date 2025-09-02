function ev = initialize_aekf_maintenance(mp, ev, jacStruct)

% Check if sim mode to avoid calling tb obj in sim mode
if mp.flagSim
    sbp_texp = mp.detector.tExpUnprobedVec; % exposure times for non-pairwise-probe images in each subband.
    psf_peaks = mp.detector.peakFluxVec;
else
    sbp_texp  = mp.tb.info.sbp_texp;
    psf_peaks = mp.tb.info.PSFpeaks;
end

% Find values to convert images back to counts rather than normalized
% intensity
ev.peak_psf_counts = zeros(1,mp.Nsbp);
ev.e_scaling = zeros(1,mp.Nsbp);

for iSubband = 1:mp.Nsbp
    % potentially set mp.detector.peakFluxVec(si) * mp.detector.tExpUnprobedVec(si) set to mp.tb.info.sbp_texp(si)*mp.tb.info.PSFpeaks(si);
    % to have cleaner setup
    ev.peak_psf_counts(iSubband) = sbp_texp(iSubband)*psf_peaks(iSubband);
    ev.e_scaling(iSubband) = sqrt(psf_peaks(iSubband));
end

% Rearrange jacobians for AEKF (3-element state)
ev.G_tot_cont = rearrange_jacobians_incoherent(mp,jacStruct,mp.dm_ind);
ev.G_tot_drift = rearrange_jacobians_incoherent(mp,jacStruct,mp.dm_drift_ind);

% Initialize EKF matrices for incoherent estimation
ev = initialize_ekf_matrices_incoherent(mp, ev, sbp_texp);

% Initialize pinned actuator check
ev.dm1.initial_pinned_actuators = mp.dm1.pinned;
if any(mp.dm_ind == 2); ev.dm2.initial_pinned_actuators = mp.dm2.pinned; end
ev.dm1.new_pinned_actuators = [];
ev.dm2.new_pinned_actuators = [];
ev.dm1.act_ele_pinned = [];
ev.dm2.act_ele_pinned = [];

% Final debug output
fprintf('=== Final AEKF Initialization State ===\n');
fprintf('x_hat size: [%d, %d]\n', size(ev.x_hat));
fprintf('P size: [%d, %d, %d, %d]\n', size(ev.P));
fprintf('Q size: [%d, %d, %d, %d]\n', size(ev.Q));
fprintf('H size: [%d, %d, %d]\n', size(ev.H));
fprintf('R size: [%d, %d, %d]\n', size(ev.R));
fprintf('G_tot_cont size: [%d, %d, %d]\n', size(ev.G_tot_cont));
fprintf('G_tot_drift size: [%d, %d, %d]\n', size(ev.G_tot_drift));
fprintf('=== End AEKF Initialization Debug ===\n');

%% =========================================================================
function G_tot = rearrange_jacobians_incoherent(mp,jacStruct,dm_inds)
% Modified to handle 3-element state per pixel

G1 = zeros(3*size(jacStruct.G1,1),mp.dm1.Nele,mp.Nsbp);
G2 = zeros(3*size(jacStruct.G2,1),mp.dm2.Nele,mp.Nsbp);

% Set up jacobian so real, imag, and incoherent components are interleaved
% and jacobian from each DM is stacked

for iSubband = 1:mp.Nsbp
    
    if any(dm_inds == 1) 
        G1_comp = jacStruct.G1(:,:,iSubband);
        G1_split = zeros(3*size(jacStruct.G1,1),mp.dm1.Nele);
        G1_split(1:3:end,:) = real(G1_comp);      % Real part
        G1_split(2:3:end,:) = imag(G1_comp);      % Imaginary part
        G1_split(3:3:end,:) = zeros(size(real(G1_comp))); % Incoherent (no DM coupling)

        G1(:,:,iSubband) = G1_split;
    else
        G1 = [];
    end

    if any(dm_inds == 2)
        G2_comp = jacStruct.G2(:,:,iSubband);
        G2_split = zeros(3*size(jacStruct.G2,1),mp.dm2.Nele);
        G2_split(1:3:end,:) = real(G2_comp);      % Real part
        G2_split(2:3:end,:) = imag(G2_comp);      % Imaginary part
        G2_split(3:3:end,:) = zeros(size(real(G2_comp))); % Incoherent (no DM coupling)

        G2(:,:,iSubband) = G2_split;
    else
        G2 = [];
    end
end

G_tot = [G1, G2];

%% =========================================================================
function ev = initialize_ekf_matrices_incoherent(mp, ev, sbp_texp)
% Modified to handle 3-element state per pixel

% Below are the definitions of the EKF matrices. There are multiple EKFs 
% defined in parallel.

%To get an idea of what the code below does, it's easier to play with 
%the toy example at https://github.com/leonidprinceton/DHMaintenanceExample

% MODIFIED: Change state size from 2 to 3 for incoherent estimation
ev.SS = 3; % Pixel state size: real(E), imag(E), incoherent_intensity
ev.BS = ev.SS*1; % EKF block size - number of pixels per EKF (currently 1). Computation time grows as the cube of BS.

ev.SL = ev.SS*mp.Fend.corr.Npix; % Total length of the state vector (all pixels).

% Debug output for initialization
fprintf('=== AEKF Initialization Debug ===\n');
fprintf('ev.SS = %d (should be 3)\n', ev.SS);
fprintf('ev.BS = %d (should be 3)\n', ev.BS);
fprintf('ev.SL = %d\n', ev.SL);
fprintf('mp.Fend.corr.Npix = %d\n', mp.Fend.corr.Npix);

%3D matrices that include all the 2D EKF matrices for all pixels at once
ev.H = zeros(floor(ev.BS/ev.SS),ev.BS,floor(ev.SL/ev.BS));
ev.R = zeros(floor(ev.BS/ev.SS),floor(ev.BS/ev.SS),floor(ev.SL/ev.BS));

%% Assemble H and R indices - MODIFIED FOR 3-ELEMENT STATE
% CORRECTED: Fix H_indices calculation for 3-element state
ev.H_indices = find(kron(eye(floor(ev.BS/ev.SS)),ones(1,ev.SS)).*ones(floor(ev.BS/ev.SS),floor(ev.BS/ev.SS)*ev.SS,floor(ev.SL/ev.BS)));
ev.R_indices = logical(eye(floor(ev.BS/ev.SS)).*ones(floor(ev.BS/ev.SS),floor(ev.BS/ev.SS),floor(ev.SL/ev.BS)));

% Debug the indices
fprintf('H_indices length: %d\n', length(ev.H_indices));
fprintf('R_indices size: [%d, %d, %d]\n', size(ev.R_indices));

% The drift covariance matrix for each pixel (or block of pixels). Needs 
% to be estimated if the model is not perfectly known.  This is a 4D
% matrix (re | im | incoherent | px | wavelength).
% Need to convert jacobian from contrast units to counts.

ev.Q = zeros(ev.SS,ev.SS,floor(ev.SL/ev.BS),mp.Nsbp);
for iSubband = 1:1:mp.Nsbp
    disp(['assembling Q for subband ',num2str(iSubband)])
   
    G_reordered = ev.G_tot_drift(:,:,iSubband);
    dm_drift_covariance = eye(size(G_reordered,2))*(mp.drift.presumed_dm_std^2);

    for i = 0:1:floor(ev.SL/ev.BS)-1
        % CORRECTED: Build Q matrix properly for 3-element state
        % Standard process noise for coherent components (first 2 elements)
        coherent_size = ev.BS-1; % Real and imaginary parts only (exclude incoherent)
        Q_coherent = G_reordered((i)*ev.BS+1:(i)*ev.BS+coherent_size,:)*dm_drift_covariance*G_reordered((i)*ev.BS+1:(i)*ev.BS+coherent_size,:).'*sbp_texp(iSubband)*(ev.e_scaling(iSubband)^2);
        
        % Assemble 3x3 Q matrix for each pixel
        Q_pixel = zeros(ev.SS,ev.SS);
        Q_pixel(1:2,1:2) = Q_coherent; % Real and imaginary components
        Q_pixel(3,3) = 1e-6; % Process noise for incoherent component (tune this value)
        
        ev.Q(:,:,i+1,iSubband) = Q_pixel;
    end
end

%% 
ev.P = ev.Q*0.0;

ev.x_hat = zeros(ev.SL,mp.Nsbp) ; % Units: sqrt(counts)
ev.x_hat0 = zeros(ev.SL,mp.Nsbp)  ;% Units: sqrt(counts)

% Initialize EKF - MODIFIED FOR 3-ELEMENT STATE
for iSubband = 1:1:mp.Nsbp
    
    % Seems to work
    % Probably doesn't make sense using probing to estimate the electric 
    % field when the photon count is low. EKF with zero initial guess 
    % should do better.
    % If an estimate is available from stroke minimization in a brighter setting:
    try %paths.E_estimated_filenames(k)
        % TODO: need to load this in main file from saved data
        E_hat = mp.est.Eest(:,iSubband) * ev.e_scaling(iSubband) * sqrt(sbp_texp(iSubband)); % assuming it is scaled (units: contrast)
    catch
        E_hat = zeros(ev.SL/ev.BS,1);%,mp.Nsbp);
    end

    % Save initial ev state - MODIFIED FOR 3-ELEMENT STATE:
    ev.x_hat0(1:3:end,iSubband) = real(E_hat);      % Real part
    ev.x_hat0(2:3:end,iSubband) = imag(E_hat);      % Imaginary part
    ev.x_hat0(3:3:end,iSubband) = zeros(size(real(E_hat))); % Initialize incoherent to zero

    % The EKF state is scaled such that the intensity is measured in photons:
    ev.x_hat(1:3:end,iSubband) = real(E_hat);       % Real part
    ev.x_hat(2:3:end,iSubband) = imag(E_hat);       % Imaginary part  
    ev.x_hat(3:3:end,iSubband) = zeros(size(real(E_hat))); % Initialize incoherent to zero
end