% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get a summed image from the back end of a single-mode optical
% fiber(s).
% 
% -------------------------------------------------------------------------
% 
% INPUTS:
% - mp: structure of model parameters.
% 
% OUTPUTS:
% - IFiberTotal: Total intensity across the bandpass from all fibers.
% 
% REVISION HISTORY:
% - Created on 4/11/2019 by Carl Coker.
%--------------------------------------------------------------------------

function IfiberTotal = falco_get_summed_image_fiber(mp)

    %--Compute the DM surfaces outside the full model to save some time
    if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
    if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
    if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end

    IfiberTotal = 0; % Initialize image

    for si=1:mp.Nsbp
        IfiberTotal = IfiberTotal +  mp.sbp_weights(si)*falco_get_sbp_image_fiber(mp,si);
    end

end %--END OF FUNCTION