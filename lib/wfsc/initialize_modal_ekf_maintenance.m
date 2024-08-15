function ev = initialize_modal_ekf_maintenance(mp, ev, jacStruct)

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

[ev.M, ev.Q] = get_DM_modes_and_drift(ev.G_tot_cont, ev.r, drift.magnitude);

ev.num_pix = size(ev.M, 2)/2;
ev.r = size(ev.M, 2);


ev.G_tot_drift = rearrange_jacobians(mp,jacStruct,mp.dm_drift_ind);

% Initialize EKF matrices
ev = initialize_ekf_matrices(mp, ev, sbp_texp);

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

function ev = initialize_ekf_matrices(mp, ev, sbp_texp)

% Below are the defnitions of the EKF matrices. There are multiple EKFs 
% defined in parallel.

 
% %% Assemble Q matrix
% ev.H_indices = find(kron(eye(floor(ev.BS/ev.SS)),ones(1,ev.SS)).*ones(floor(ev.BS/ev.SS),floor(ev.BS/ev.SS)*ev.SS,floor(ev.SL/ev.BS)));
% ev.R_indices = logical(eye(floor(ev.BS/ev.SS)).*ones(floor(ev.BS/ev.SS),floor(ev.BS/ev.SS),floor(ev.SL/ev.BS)));
% 
% % The drift covariance matrix for each pixel (or block of pixels). Needs 
% % to be estimated if the model is not perfectly known.  This is a 4D
% % matrix (re | im | px | wavelength).
% % Need to convert jacobian from contrast units to counts.
% 
ev.R = zeros(ev.num_pix, ev.num_pix);
ev.H = zeros(ev.num_pix, ev.r);

ev.E_floor = zeros(ev.r, mp.Nsbp);

for iSubband = 1:1:mp.Nsbp
    disp(['assembling Q for subband ',num2str(iSubband)])
   
    G_reordered = ev.G_tot_drift(:,:,iSubband);
    dm_drift_covariance = eye(size(G_reordered,2))*(mp.drift.presumed_dm_std^2);

    ev.Q(:,:,iSubband) = G_reordered*dm_drift_covariance*G_reordered.'*sbp_texp(iSubband)*(ev.e_scaling(iSubband)^2);
end

%% 
ev.P = ev.Q*0.0;

ev.x_hat = zeros(ev.r,mp.Nsbp) ; % Units: sqrt(counts)
ev.x_hat0 = zeros(ev.r,mp.Nsbp)  ;% Units: sqrt(counts)

% Initialized EKF
for iSubband = 1:1:mp.Nsbp
    
    % Seems to work
    % Probably doesn't make sense using probing to estimate the electric 
    % field when the photon count is low. EKF with zero initial guess 
    % should do better.
    % If an estimate is availible from stroke minimization in a brighter setting:
    try %paths.E_estimated_filenames(k)
        % TODO: need to load this in main file from saved data
        % E_hat = mp.est.Eest(:,iSubband) * ev.e_scaling(iSubband) * sqrt(mp.tb.info.sbp_texp(iSubband)); % assuming it is scaled (units: contrast)
        E_hat = ev.E_floor + mp.est.Eest(:, iSubband) * ev.e_scaling(iSubband) * sqrt(mp.tb.info.sbp_texp(iSubband));
    catch
        E_hat = zeros(ev.r,1);%,mp.Nsbp);
    end

    % Save initial ev state:
    ev.x_hat0(1:ev.r,iSubband) = E_hat * ev.e_scaling(iSubband) * sqrt(sbp_texp(iSubband));

    % The EKF state is scaled such that the intensity is measured in photons:
    ev.x_hat(1:ev.r,iSubband) = ev.xhat0(1:ev.r, iSubband);
 end

end

function [M, QM] = get_DM_modes_and_drift(G, r, drift)
    [U, S, V] = svd(G, 'econ');
    
    S = diag(S);
    if isempty(r)
        r = length(diag(S));
    end

    M = V(:, 1:r);

    QM = diag(S(1:r))*(U(:, 1:r)')*(drift^2.*eye(size(G, 1))) * U(:, 1:r) * diag(S(1:r));
end
