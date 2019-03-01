% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Isum = falco_get_summed_image(mp)
%
% Function to get a broadband image over the entire bandpass by summing the 
% sub-bandpass images.
%
%--INPUTS
% mp = structure of all model parameters
%
%--OUTPUTS
% Ibandavg = band-averaged image in units of normalized intensity
%
%--REVISION HISTORY
% - Simplified on 2019-03-01 by A.J. Riggs to loop over falco_get_sbp_image.m 
%--------------------------------------------------------------------------

function Ibandavg = falco_get_summed_image(mp)

    %--Compute the DM surfaces outside the full model to save some time
    if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
    if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
    if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end

    % facContrastToCounts = model_params.texp*model_params.peakCountsPerPixPerSec;
    Ibandavg = 0; % Initialize image

    for si=1:mp.Nsbp    
        Ibandavg = Ibandavg +  mp.sbp_weights(si)*falco_get_sbp_image(mp,si);
    end

end %--END OF FUNCTION

