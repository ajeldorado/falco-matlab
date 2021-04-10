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


    IfiberTotal = 0; % Initialize image

    for si=1:mp.Nsbp
        Ifiber = falco_get_sbp_image_fiber(mp,si);
        IfiberTotal = IfiberTotal +  mp.sbp_weights(si)*Ifiber;
    end

end %--END OF FUNCTION