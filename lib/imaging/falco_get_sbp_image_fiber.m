% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass from an optical
% fiber.
% 
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - si = index of sub-bandpass for which to take the image
%
% OUTPUTS
% - Im: sub-bandpass image in units of normalized intensity
%
% REVISION HISTORY
% - Created on 4/11/2019 by Carl Coker.

function ImNI = falco_get_sbp_image_fiber(mp,si)

    if(mp.flagSim) %--Generate simulated image
        ImNI = falco_get_sim_sbp_image_fiber(mp,si);
    else
        error('Testbed functionality not implemented for fibers yet.')
    end

end %--END OF FUNCTION