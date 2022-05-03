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
% Ecube = Cube of the 2-D stellar complex field at each wavelength

function Ecube = falco_get_Efield_cube(mp)

    for si = mp.Nsbp:-1:1 
        modvar = ModelVariables;
        modvar.sbpIndex   = si;
        modvar.wpsbpIndex = mp.wi_ref;
        modvar.whichSource = 'star';
        modvar.starIndex = 1;
        Ecube(:,:,si) = model_full(mp, modvar);
    end

end
