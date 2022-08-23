% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Flesh out the rest of the workspace prior to wavefront estimation and control.
%
% [mp, out] = falco_flesh_out_workspace(mp)

function [mp, out] = falco_flesh_out_workspace(mp)

mp = falco_set_optional_variables(mp);
mp = falco_verify_key_values(mp);

mp = falco_set_spectral_properties(mp);
mp = falco_set_jacobian_modal_weights(mp);

%--Pupil Masks
mp = falco_compute_entrance_pupil_coordinates(mp);
mp = falco_compute_apodizer_shape(mp);
mp = falco_crop_lyot_stop(mp);
mp = falco_compute_lyot_stop_coordinates(mp);
falco_plot_superimposed_pupil_masks(mp); %--Visually inspect relative pupil mask alignment

%--Focal plane masks
mp = falco_gen_fpm(mp);
mp = falco_compute_fpm_coordinates(mp); 

%--Final focal plane
mp = falco_compute_Fend_resolution(mp); 
mp = falco_configure_dark_hole_region(mp);
mp = falco_set_jacobian_spatial_weights(mp);

%--Wavefront sensor
mp = falco_configure_wfs(mp);

%--DM1 and DM2
mp = falco_configure_dm1_and_dm2(mp);
mp = falco_gen_dm_stops(mp);
mp = falco_set_dm_surface_padding(mp);

mp = falco_set_initial_Efields(mp);
mp = falco_compute_psf_norm_factor(mp);

out = falco_init_storage_arrays(mp); % Initialize Arrays to Store Performance History

%--Save the config file
fn_config = [mp.path.config filesep mp.runLabel,'_config.mat'];
save(fn_config, 'mp')
fprintf('Saved the config file: \t%s\n', fn_config)

end %--END OF FUNCTION
