% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% summedImage = falco_get_summed_image(mp)
%
% Get a broadband image over the entire bandpass by summing subband images.
%
% INPUTS
% ------
% mp : structure of all model parameters
%
% OUTPUTS
% -------
% summedImage : PSF averaged over all subbands and polarization states [normalized intensity]

function summedImage = falco_get_summed_image(mp)

    summedImage = 0;
    for iSubband = 1:mp.Nsbp    
        summedImage = summedImage +  mp.sbp_weights(iSubband)*falco_get_sbp_image(mp, iSubband);
    end

end
