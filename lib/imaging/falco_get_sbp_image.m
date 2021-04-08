% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Take an image in the specified subband.
%
% INPUTS
% ------
% mp : structure of model parameters
% iSubband : index of subband for which to take the image
%
% OUTPUTS
% -------
% subbandImage : subband image in units of normalized intensity

function subbandImage = falco_get_sbp_image(mp, iSubband)

    if mp.flagSim
        subbandImage = falco_get_sim_sbp_image(mp, iSubband);
    else
        subbandImage = falco_get_testbed_sbp_image(mp, iSubband);
    end

end
