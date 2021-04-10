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
        %--Compute the DM surfaces outside the full model to save some time
        if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
        if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
        if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end
        ImNI = falco_get_sim_sbp_image_fiber(mp,si);
        ImNI = sum(sum(ImNI));
    else
        ImNI = falco_get_hcst_sbp_image_fiber(mp,si);
    end

end %--END OF FUNCTION
