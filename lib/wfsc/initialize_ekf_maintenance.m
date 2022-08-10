function ev = initialize_ekf_maintenance(mp, ev, jacStruct)


% Find values to convert images back to counts rather than normalized
% intensity
ev.peak_psf_counts = zeros(1,mp.Nsbp);
ev.e_scaling = zeros(1,mp.Nsbp);

%     sbp_texp  = tb.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
%     numCoadds = tb.info.sbp_numCoadds(si);% Number of coadds to use for each sub-bandpass
%     PSFpeak   = tb.info.PSFpeaks(si);% counts per second 
%     PSFpeak_counts = PSFpeak*sbp_texp*numCoadds;

for iSubband = 1:mp.Nsbp

    % potentially set mp.detector.peakFluxVec(si) * mp.detector.tExpUnprobedVec(si) set to mp.tb.info.sbp_texp(si)*mp.tb.info.PSFpeaks(si);
    % to have cleaner setup
    ev.peak_psf_counts(iSubband) = mp.tb.info.sbp_texp(iSubband)*mp.tb.info.PSFpeaks(iSubband);
    ev.e_scaling(iSubband) = sqrt(mp.tb.info.PSFpeaks(iSubband));

end

% Get e_scaling
% TODO: need e_scaling

% Rearrange jacobians
ev = rearrange_jacobians(mp,ev,jacStruct);


% Initialize EKF matrices
ev = initialize_ekf_matrices(mp, ev);


% Initialize pinned actuator check
ev.dm1.new_pinned_actuators = [];
ev.dm2.new_pinned_actuators = [];


end


function ev = rearrange_jacobians(mp,ev,jacStruct)

% active_dms = ismember([1:9],mp.dm_ind);

% assume we are max using 2 DMs for estimation
% jacStruct.G_tot = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele*active_dms(1) + mp.dm2.Nele*active_dms(2),Nsbp);

G1 = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele,mp.Nsbp);
G2 = zeros(2*size(jacStruct.G1,1),mp.dm2.Nele,mp.Nsbp);

% Set up jacobian so real and imag components alternate and jacobian from
% each DM is stacked
for iSubband = 1:mp.Nsbp
    
    if any(mp.dm_ind == 1)
        G1_comp = jacStruct.G1(:,:,iSubband);
        G1_split = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele);
        G1_split(1:2:end,:) = real(G1_comp);
        G1_split(2:2:end,:) = imag(G1_comp);

        G1(:,:,iSubband) = G1_split;

    else
        G1 = [];
    end

    if any(mp.dm_ind == 2)
        G2_comp = jacStruct.G2(:,:,iSubband);
        G2_split = zeros(2*size(jacStruct.G2,1),mp.dm2.Nele);
        G2_split(1:2:end,:) = real(G2_comp);
        G2_split(2:2:end,:) = imag(G2_comp);

        G2(:,:,iSubband) = G2_split;

    else
        G2 = [];

    end
end

ev.G_tot = [G1, G2];


end

function ev = initialize_ekf_matrices(mp, ev)

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

% % TODO: check if we have 2 DMs available, consider injecting drift only on
% % one
ev.Q = zeros(ev.SS,ev.SS,floor(ev.SL/ev.BS),mp.Nsbp);
for iSubband = 1:1:mp.Nsbp
    disp(['assembling Q for subband ',num2str(iSubband)])
   
    G_reordered = ev.G_tot(:,:,iSubband);
    dm_drift_covariance = eye(size(G_reordered,2))*(mp.drift.presumed_dm_std^2);

    for i = 0:1:floor(ev.SL/ev.BS)-1
%         size(dm_drift_covariance)
%         size(G_reordered((i)*ev.BS+1:(i+1)*ev.BS,:))
        ev.Q(:,:,i+1,iSubband) = G_reordered((i)*ev.BS+1:(i+1)*ev.BS,:)*dm_drift_covariance*G_reordered(i*ev.BS+1:(i+1)*ev.BS,:).'*mp.tb.info.sbp_texp(iSubband)*(ev.e_scaling(iSubband)^2);
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
        E_hat = mp.est.Eest(:,iSubband) * ev.e_scaling(iSubband) * sqrt(mp.tb.info.sbp_texp(iSubband)); % assuming it is scaled (units: contrast)
    catch
        E_hat = zeros(ev.SL/ev.BS,1);%,mp.Nsbp);
    end


    % UPDATE EXPOSURE TIME HERE ? ****************************************
    % Save initial ev state:
    ev.x_hat0(1:ev.SS:end,iSubband) = real(E_hat) * ev.e_scaling(iSubband) * sqrt(mp.tb.info.sbp_texp(iSubband));
    ev.x_hat0(2:ev.SS:end,iSubband) = imag(E_hat) * ev.e_scaling(iSubband) * sqrt(mp.tb.info.sbp_texp(iSubband));

    % Optional: save initial state for debugging
    % hicat.util.write_fits(x_hat0, os.path.join(initial_path, f"x_hat0_{int(wavelength)}nm.fits"))

    % The EKF state is scaled such that the intensity is measured in photons:
    ev.x_hat(1:ev.SS:end,iSubband) = real(E_hat) * ev.e_scaling(iSubband) * sqrt(mp.tb.info.sbp_texp(iSubband));
    ev.x_hat(2:ev.SS:end,iSubband) = imag(E_hat) * ev.e_scaling(iSubband) * sqrt(mp.tb.info.sbp_texp(iSubband));
 end

end



% TODO: these probable don't do what I want them to do right now

% function subbandImage = falco_get_sim_sbp_peak_counts(mp, iSubband)
% 
%     %--Compute the DM surfaces outside the full model to save time
%     if any(mp.dm_ind == 1); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1, mp.dm1.dx, mp.dm1.NdmPad); end
%     if any(mp.dm_ind == 2); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2, mp.dm2.dx, mp.dm2.NdmPad); end
%     if any(mp.dm_ind == 9); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9, mp.dm9); end
% 
%     %--Loop over all wavelengths and polarizations
%     Npol = length(mp.full.pol_conds);
%     indexComboArray = allcomb(1:mp.Nwpsbp, 1:Npol, 1:mp.star.count).'; % dimensions: [3 x mp.Nwpsbp*Npol*mp.star.count ]
%     Ncombos = size(indexComboArray, 2);
% 
%     if mp.flagParfor
%         parfor iCombo = 1:Ncombos
%             Iall{iCombo} = falco_compute_subband_image_peak(mp, indexComboArray, iCombo, iSubband);
%         end
%     else
%         for iCombo = Ncombos:-1:1
%             Iall{iCombo} = falco_compute_subband_image_peak(mp, indexComboArray, iCombo, iSubband);
%         end
%     end
% 
%     %--Apply the spectral weights and sum
%     subbandImage = 0; 
%     for iCombo = 1:Ncombos
%         subbandImage = subbandImage + Iall{iCombo};  
%     end
%     
%     if mp.flagImageNoise
%         subbandImage = falco_add_noise_to_subband_image(mp, subbandImage, iSubband);
%     end
% 
% 
% end %--END OF FUNCTION
% 
% 
% function Iout = falco_compute_subband_image_peak(mp, indexComboArray, iCombo, iSubband)
% 
%     % Generate the weighted, normalized intensity image for a single
%     % wavelength, polarization, and star in the specified subband.
% 
%     % Extract indices
%     iWavelength = indexComboArray(1, iCombo);
%     iPol = indexComboArray(2, iCombo);
%     iStar = indexComboArray(3, iCombo);
% 
%     % Compute E-field
%     modvar = ModelVariables;
%     modvar.whichSource = 'star';
%     modvar.sbpIndex = iSubband;
%     modvar.wpsbpIndex = iWavelength;
%     modvar.starIndex = iStar;
%     mp.full.polaxis = mp.full.pol_conds(iPol); % used only in PROPER full models
%     Estar = model_full(mp, modvar);
% 
%     % Apply wavelength weight within subband.
%     % Assume polarizations are evenly weighted.
%     Iout = length(mp.full.pol_conds) * abs(Estar).^2;
%     
% end %--END OF FUNCTION
