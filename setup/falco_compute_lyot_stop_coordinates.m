% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute the resolution and coordinates for the cropped-down Lyot stop.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_compute_lyot_stop_coordinates(mp)

    % Full model
    if ~mp.full.flagPROPER
        mp.P4.full.dx = mp.P4.D/mp.P4.full.Nbeam; % [meters per pixel]
    end

    % Compact model
    mp.P4.compact.dx = mp.P4.D/mp.P4.compact.Nbeam; % [meters per pixel]
    if strcmpi(mp.centering, 'pixel')
        mp.P4.compact.xs = (-mp.P4.compact.Narr/2:(mp.P4.compact.Narr/2-1))*mp.P4.compact.dx;
    elseif strcmpi(mp.centering, 'interpixel')
        mp.P4.compact.xs = (-(mp.P4.compact.Narr-1)/2:(mp.P4.compact.Narr-1)/2)*mp.P4.compact.dx;
    end
    mp.P4.compact.ys = mp.P4.compact.xs.';

end
