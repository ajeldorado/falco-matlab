% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass.
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
% - Created on 2019-02-06 by A.J. Riggs.

function ImNI = falco_get_sbp_image(mp,si)

    if(mp.flagSim) %--Generate simulated image
        ImNI = falco_get_sim_sbp_image(mp,si);
    else %--Retrieve testbed image
        ImNI = falco_get_testbed_sbp_image(mp,si);
    end

end %--END OF FUNCTION