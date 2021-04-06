% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute the array size at the apodizer.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_compute_apodizer_coordinates(mp)

    if isfield(mp.P3.compact, 'mask')
        mp.P3.compact.mask = pad_to_even_square(mp.P3.compact.mask); % force to be square array
        mp.P3.compact.Narr = length(mp.P3.compact.mask);
    end
    
    if isfield(mp.P3.full, 'mask')
        mp.P3.full.mask = pad_to_even_square(mp.P3.full.mask); % force to be square array
        mp.P3.full.Narr = length(mp.P3.full.mask);
    end

end
