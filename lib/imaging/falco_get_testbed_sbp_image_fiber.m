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
% - normI: Normalized intensity in the sub-bandpass
%          	(i.e. approximate raw contrast but normalized 
%           by a photometry measurement at a single offset)
%
% REVISION HISTORY
% - Created on 2020-02-11 by Carl Coker.

function normI = falco_get_testbed_sbp_image_fiber(mp,si)

    switch upper(mp.testbed)
        case 'HCST'
            normI = falco_get_hcst_sbp_image_fiber(mp,si);
        %Can put other testbeds or fancy models here
        %case 'WHATEVER'
            %normI = falco_get_whatever_sbp_image_fiber(mp,si); 
        otherwise
            error('Case not recognized.  Check falco_get_testbed_sbp_image.m to make sure this case is specified.');
    end

end %--END OF FUNCTION