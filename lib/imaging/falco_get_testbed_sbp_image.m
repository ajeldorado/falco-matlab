% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass from a testbed. 
% This function calls an equivalent sub-function depending on mp.testbed. 
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
% - Created on 2019-03-22 by G. Ruane
% - Silent failure if no testbed case sepcified removed by C. Coker
% (2019-04-11)

function normI = falco_get_testbed_sbp_image(mp,si)
    
    switch upper(mp.testbed)
        case 'HCST'
            normI = falco_get_hcst_sbp_image(mp,si);
        case 'GPCT'
            normI = falco_get_gpct_sbp_image(mp,si); 
        %Can put other testbeds or fancy models here
        %case 'WHATEVER'
            %normI = falco_get_whatever_sbp_image(mp,si); 
        otherwise
            error('Case not recognized.  Check falco_get_testbed_sbp_image.m to make sure this case is specified.');
    end
    
end