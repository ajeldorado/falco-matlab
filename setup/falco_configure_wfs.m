% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Configure the Zernike wavefront sensor (ZWFS).
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_configure_wfs(mp)

    if mp.flagWFS
        
        mp = falco_gen_fpm_zwfs(mp);
        mp = falco_compute_fpm_zwfs_coordinates(mp);
        mp = falco_set_wfs_spectral_properties(mp);
        
        % Intialize E-field seen by WFS (will be different to coronagraph 
        % if using WFS out of band). This is updated in model_ZWFS.
        mp.wfs.E = ones(mp.P1.full.Narr, mp.P1.full.Narr, mp.wfs.Nwpsbp, mp.wfs.Nsbp);
        
    end
        
end %--END OF FUNCTION