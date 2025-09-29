function ev = initialize_ekf_maintenance(mp, ev, jacStruct)

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



% Rearrange jacobians
ev.G_tot_cont = rearrange_jacobians(mp,jacStruct,mp.dm_ind);
ev.G_tot_drift = rearrange_jacobians(mp,jacStruct,mp.dm_drift_ind);

% Initialize EKF matrices
switch lower(mp.estimator)
    case{'ekf_maintenance'}
        ev = initialize_ekf_matrices(mp, ev, sbp_texp);
    case{'modal_ekf_maintenance'}
        ev = initialize_modal_ekf_matrices(mp, ev);
end

% Initialize pinned actuator check
ev.dm1.initial_pinned_actuators = mp.dm1.pinned;
if any(mp.dm_ind == 2); ev.dm2.initial_pinned_actuators = mp.dm2.pinned; end
ev.dm1.new_pinned_actuators = [];
ev.dm2.new_pinned_actuators = [];
ev.dm1.act_ele_pinned = [];
ev.dm2.act_ele_pinned = [];

end


function G_tot = rearrange_jacobians(mp,jacStruct,dm_inds)

G1 = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele,mp.Nsbp);
G2 = zeros(2*size(jacStruct.G2,1),mp.dm2.Nele,mp.Nsbp);

% Set up jacobian so real and imag components alternate and jacobian from
% each DM is stacked

for iSubband = 1:mp.Nsbp
    
    if any(dm_inds == 1) 
        G1_comp = jacStruct.G1(:,:,iSubband);
        G1_split = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele);
        G1_split(1:2:end,:) = real(G1_comp);
        G1_split(2:2:end,:) = imag(G1_comp);

        G1(:,:,iSubband) = G1_split;

    else
        G1 = [];
    end

    if any(dm_inds == 2)
        G2_comp = jacStruct.G2(:,:,iSubband);
        G2_split = zeros(2*size(jacStruct.G2,1),mp.dm2.Nele);
        G2_split(1:2:end,:) = real(G2_comp);
        G2_split(2:2:end,:) = imag(G2_comp);

        G2(:,:,iSubband) = G2_split;

    else
        G2 = [];

    end
end

G_tot = [G1, G2];


end

function ev = initialize_modal_ekf_matrices(mp, ev)

    fprintf('Initializing modal EKF...\n');
    if isfield(mp.est, 'r')
        Nact_est = mp.est.r;
        Nact_drift = mp.est.r;
    else
        nele_vals = [mp.dm1.Nele, mp.dm2.Nele];
        Nact_est = sum(nele_vals(mp.dm_ind));
        Nact_drift = sum(nele_vals(mp.dm_drift_ind));
    end

    %--Initial DM command, based on all available DM actuators
    ev.SS = 2; % Pixel state size. Two for real and imaginary parts of the electric field. If incoherent intensity is not ignored, SS should be 3 and the EKF modified accordingly.
    ev.x_hat = zeros(Nact_est, mp.Nsbp);
    ev.Q = zeros(Nact_est,Nact_est,mp.Nsbp);
    ev.P = zeros(Nact_est,Nact_est,mp.Nsbp);
    for iSubband = 1:1:mp.Nsbp
        %--Initial covariance, based on drifting DM actuators
        act_drift = 1:Nact_drift;
        ev.Q(act_drift,act_drift,iSubband) = eye(Nact_drift) * mp.drift.presumed_dm_std^2;
        
        %--Process noise covariance, based on drifting DM actuators
        ev.P(act_drift,act_drift,iSubband) = eye(Nact_drift) * mp.drift.presumed_dm_std^2; %+ eye(size(mp.dm_drift_ind, 2) * Nact) * mp.est.testbed_drift^2;
    end

end

function ev = initialize_ekf_matrices(mp, ev, sbp_texp)

% Below are the defnitions of the EKF matrices. There are multiple EKFs 
% defined in parallel.

%To get an idea of what the code below does, it's easier to play with 
%the toy example at https://github.com/leonidprinceton/DHMaintenanceExample
ev.SS = 2; % Pixel state size. Two for real and imaginary parts of the electric field. If incoherent intensity is not ignored, SS should be 3 and the EKF modified accordingly.
ev.BS = ev.SS*1; % EKF block size - number of pixels per EKF (currently 1). Computation time grows as the cube of BS.


ev.SL = ev.SS*mp.Fend.corr.Npix;%sum(size(jacStruct.G1,1));%Total length of the sate vector (all pixels).

%3D matrices that include all the 2D EKF matrices for all pixels at once
ev.H = zeros(floor(ev.BS/ev.SS),ev.BS,floor(ev.SL/ev.BS));

ev.R = zeros(floor(ev.BS/ev.SS),floor(ev.BS/ev.SS),floor(ev.SL/ev.BS));

%% Assemble Q matrix
ev.H_indices = find(kron(eye(floor(ev.BS/ev.SS)),ones(1,ev.SS)).*ones(floor(ev.BS/ev.SS),floor(ev.BS/ev.SS)*ev.SS,floor(ev.SL/ev.BS)));
ev.R_indices = logical(eye(floor(ev.BS/ev.SS)).*ones(floor(ev.BS/ev.SS),floor(ev.BS/ev.SS),floor(ev.SL/ev.BS)));



% The drift covariance matrix for each pixel (or block of pixels). Needs 
% to be estimated if the model is not perfectly known.  This is a 4D
% matrix (re | im | px | wavelength).
% Need to convert jacobian from contrast units to counts.

ev.Q = zeros(ev.SS,ev.SS,floor(ev.SL/ev.BS),mp.Nsbp);
for iSubband = 1:1:mp.Nsbp
    disp(['assembling Q for subband ',num2str(iSubband)])
   
    G_reordered = ev.G_tot_drift(:,:,iSubband);
    dm_drift_covariance = eye(size(G_reordered,2))*(mp.drift.presumed_dm_std^2);

    for i = 0:1:floor(ev.SL/ev.BS)-1
        ev.Q(:,:,i+1,iSubband) = G_reordered((i)*ev.BS+1:(i+1)*ev.BS,:)*dm_drift_covariance*G_reordered(i*ev.BS+1:(i+1)*ev.BS,:).'*sbp_texp(iSubband)*(ev.e_scaling(iSubband)^2);
    end
end

%% 
ev.P = ev.Q*0.0;

ev.x_hat = zeros(ev.SL,mp.Nsbp) ; % Units: sqrt(counts)
ev.x_hat0 = zeros(ev.SL,mp.Nsbp)  ;% Units: sqrt(counts)

% Initialized EKF
for iSubband = 1:1:mp.Nsbp
    
    % Seems to work
    % Probably doesn't make sense using probing to estimate the electric 
    % field when the photon count is low. EKF with zero initial guess 
    % should do better.
    % If an estimate is availible from stroke minimization in a brighter setting:
    try %paths.E_estimated_filenames(k)
        % TODO: need to load this in main file from saved data
        E_hat = mp.est.Eest(:,iSubband); % * ev.e_scaling(iSubband) * sqrt(mp.tb.info.sbp_texp(iSubband)); % assuming it is scaled (units: contrast)
    catch
        E_hat = zeros(ev.SL/ev.BS,1);%,mp.Nsbp);
    end

    % Save initial ev state:
    ev.x_hat0(1:ev.SS:end,iSubband) = real(E_hat) * ev.e_scaling(iSubband) * sqrt(sbp_texp(iSubband));
    ev.x_hat0(2:ev.SS:end,iSubband) = imag(E_hat) * ev.e_scaling(iSubband) * sqrt(sbp_texp(iSubband));

    % The EKF state is scaled such that the intensity is measured in photons:
    ev.x_hat(1:ev.SS:end,iSubband) = real(E_hat) * ev.e_scaling(iSubband) * sqrt(sbp_texp(iSubband));
    ev.x_hat(2:ev.SS:end,iSubband) = imag(E_hat) * ev.e_scaling(iSubband) * sqrt(sbp_texp(iSubband));
 end

end