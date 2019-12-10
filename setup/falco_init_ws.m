% Copyright 2018,2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to finish initializing the workspace prior to wavefront
% estimation and control.

function [mp,out] = falco_init_ws(fn_config)

load(fn_config,'mp'); %% Read inputs as structures from a .mat config file

mp = falco_set_optional_variables(mp); %% Optional/hidden boolean flags and variables

%--File Paths for Data Storage (excluded from git)
if(isfield(mp.path,'ws')==false); mp.path.ws = [mp.path.falco 'data' filesep 'ws' filesep]; end % Store final workspace data here
if(isfield(mp.path,'maps')==false); mp.path.falcoaps = [mp.path.falco 'maps' filesep]; end % Maps go here
if(isfield(mp.path,'jac')==false); mp.path.jac = [mp.path.falco 'data' filesep 'jac' filesep]; end % Store the control Jacobians here
if(isfield(mp.path,'images')==false); mp.path.images = [mp.path.falco 'data' filesep 'images' filesep]; end % Store all full, reduced images here
if(isfield(mp.path,'dm')==false); mp.path.dm = [mp.path.falco 'data' filesep 'DM' filesep]; end % Store DM command maps here
if(isfield(mp.path,'ws_inprogress')==false); mp.path.ws_inprogress = [mp.path.falco 'data' filesep 'ws_inprogress' filesep]; end % Store in progress workspace data here

mp.mas2lam0D = 1/(mp.lambda0/mp.P1.D*180/pi*3600*1000); %% Conversion factor: milliarcseconds (mas) to lambda0/D

mp = falco_set_spectral_weights(mp);
mp = falco_set_jacobian_weights(mp); %% Set Weighting of the Control Jacobian (of Zernike Modes and Subbands)

%--Pupil Masks
mp = falco_gen_chosen_pupil(mp);
mp = falco_gen_chosen_apodizer(mp);
mp = falco_gen_chosen_lyot_stop(mp);
falco_plot_superposed_pupil_masks(mp); %--Visually inspect relative pupil mask alignment

mp = falco_configure_dm1_to_dm9(mp); %% Initialize some basic attributes for all DMs.

mp = falco_gen_FPM(mp); %% Generate FPM if necessary
mp = falco_get_FPM_coordinates(mp); 

%%--Sampling/Resolution and Scoring/Correction Masks for Final Focal Plane (Fend)
mp.Fend.dxi = (mp.fl*mp.lambda0/mp.P4.D)/mp.Fend.res; % sampling at Fend.[meters]
mp.Fend.deta = mp.Fend.dxi; % sampling at Fend.[meters]    
if(mp.flagLenslet)
    mp.Fend.lenslet.D = 2*mp.Fend.res*mp.Fend.lensletWavRad*mp.Fend.dxi;
    mp.Fend.x_lenslet_phys = mp.Fend.dxi*mp.Fend.res*mp.Fend.x_lenslet;
    mp.Fend.y_lenslet_phys = mp.Fend.deta*mp.Fend.res*mp.Fend.y_lenslet;

    mp.F5.dxi = mp.lensletFL*mp.lambda0/mp.Fend.lenslet.D/mp.F5.res;
    mp.F5.deta = mp.F5.dxi;
end
    
mp = falco_configure_dark_hole_region(mp); %% Software Mask for Correction (corr) and Scoring (score)

%%--Spatial weighting of dark hole pixel intensity. 
if(mp.flagFiber && mp.flagLenslet)
    mp.WspatialVec = ones(mp.Fend.Nlens,1);
else
    mp = falco_set_spatial_weights(mp);
    mp.WspatialVec = mp.Wspatial(mp.Fend.corr.maskBool);    
end

mp = falco_configure_dm1_and_dm2(mp); %% Flesh out the dm1 and dm2 structures

%%--DM Aperture Masks (moved here because the commands mp.dm2.compact = mp.dm2; and mp.dm1.compact = mp.dm1; otherwise would overwrite the compact model masks)
if(mp.flagDM1stop)
    mp.dm1.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm1.Dstop,mp.centering);
    mp.dm1.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm1.Dstop,mp.centering);
end
if(mp.flagDM2stop)
    mp.dm2.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm2.Dstop,mp.centering);
    mp.dm2.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm2.Dstop,mp.centering);
end

mp = falco_set_dm_surface_padding(mp); %% DM Surface Array Sizes for Angular Spectrum Propagation with FFTs

%%--Initial Electric Fields for Star and Exoplanet
if(isfield(mp.P1.full,'E')==false)
    mp.P1.full.E  = ones(mp.P1.full.Narr,mp.P1.full.Narr,mp.Nwpsbp,mp.Nsbp);
end
if(isfield(mp,'Eplanet')==false)
    mp.Eplanet = mp.P1.full.E; %--Phase ramp added later in propagation model
end
if(isfield(mp.P1.compact,'E')==false)
    mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp);
end
mp.sumPupil = sum(sum(abs(mp.P1.compact.mask.*padOrCropEven(mean(mp.P1.compact.E,3),size(mp.P1.compact.mask,1) )).^2));

% mp.Fend.mask = ones(mp.Fend.Neta,mp.Fend.Nxi); %% Field Stop at Fend (as a software mask) (NOT INCLUDED YET)

%%--Contrast to Normalized Intensity Map Calculation (NOT INCLUDED YET)

mp = falco_get_PSF_norm_factor(mp); %% Get the starlight normalization factor for the compact and full models (to convert images to normalized intensity)

out = falco_init_storage_arrays(mp); %% Initialize Arrays to Store Performance History

fprintf('\nBeginning Trial %d of Series %d.\n',mp.TrialNum,mp.SeriesNum);

end %--END OF FUNCTION