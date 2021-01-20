% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Imean = falco_get_Efield_cube(mp)
%
% Function to get a broadband image over the entire bandpass by summing the 
% sub-bandpass images.
%
%--INPUTS
% mp = structure of all model parameters
%
%--OUTPUTS
% Ecube = 2-D stellar complex field at each wavelength
%
%--REVISION HISTORY
% - Modified from falco_get_summed_image.m on 2020-12-23 by G. Ruane 
% - Modified on 2019-05-06 by A.J. Riggs to include a different option for
% looping over a full model in PROPER.
% - Simplified on 2019-03-01 by A.J. Riggs to loop over falco_get_sbp_image.m 
%--------------------------------------------------------------------------

function Ecube = falco_get_Efield_cube(mp)

    for si=1:mp.Nsbp    
        modvar.sbpIndex   = si;
        modvar.wpsbpIndex = mp.wi_ref;
        modvar.whichSource = 'star';
        modvar.starIndex = 1;
        Ecube(:,:,si) = model_full(mp, modvar);
    end

end %--END OF FUNCTION
