% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Flesh out the rest of the workspace prior to wavefront estimation and control.

function [mp, out] = falco_flesh_out_workspace(mp)

mp = falco_set_optional_variables(mp);

mp = falco_set_spectral_properties(mp);
mp = falco_set_jacobian_modal_weights(mp);

%--Pupil Masks
mp = falco_gen_chosen_pupil(mp);
mp = falco_gen_chosen_apodizer(mp);
mp = falco_gen_chosen_lyot_stop(mp);
falco_plot_superimposed_pupil_masks(mp); %--Visually inspect relative pupil mask alignment

%--Focal planes
mp = falco_gen_FPM(mp);
mp = falco_get_FPM_coordinates(mp); 
mp = falco_get_Fend_resolution(mp); 
mp = falco_configure_dark_hole_region(mp);
mp = falco_set_jacobian_spatial_weights(mp);

%--Wavefront sensor
mp = falco_configure_WFS(mp);

%--DM1 and DM2
mp = falco_configure_dm1_and_dm2(mp);
mp = falco_gen_DM_stops(mp);
mp = falco_set_dm_surface_padding(mp);

mp = falco_set_initial_Efields(mp);

mp = falco_compute_PSF_norm_factor(mp);

out = falco_init_storage_arrays(mp); %% Initialize Arrays to Store Performance History

%--Save the config file
fn_config = [mp.path.config mp.runLabel,'_config.mat'];
save(fn_config,'mp')
fprintf('Saved the config file: \t%s\n',fn_config)

fprintf('\nBeginning Trial %d of Series %d.\n',mp.TrialNum,mp.SeriesNum);

end %--END OF FUNCTION