% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute the resolution and coordinates at the entrance pupil (plane P1).
% Values also are true at P2 in FALCO models.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_compute_entrance_pupil_coordinates(mp)

    % Compact model: Make sure mask is square
    mp.P1.compact.mask = pad_to_even_square(mp.P1.compact.mask); % force to be square

    % Compact model coordinates normalized to pupil diameter
    % Used to make the tip/tilted input wavefront within the compact model.
    mp.P1.compact.Narr = length(mp.P1.compact.mask);
    if strcmpi(mp.centering, 'pixel')
        mp.P2.compact.xsDL = (-mp.P1.compact.Narr/2:(mp.P1.compact.Narr/2-1))*mp.P2.compact.dx/mp.P2.D;
    elseif strcmpi(mp.centering, 'interpixel')
        mp.P2.compact.xsDL = (-(mp.P1.compact.Narr-1)/2:(mp.P1.compact.Narr-1)/2)*mp.P2.compact.dx/mp.P2.D;
    end
    [mp.P2.compact.XsDL, mp.P2.compact.YsDL] = meshgrid(mp.P2.compact.xsDL);

    % Full model: Number of points across array
    if mp.full.flagPROPER
        if strcmpi(mp.centering, 'pixel')
            mp.P1.full.Narr = ceil_even(mp.P1.full.Nbeam+1);
        elseif strcmpi(mp.centering, 'interpixel')
            mp.P1.full.Narr = ceil_even(mp.P1.full.Nbeam);
        end
    else
        mp.P1.full.mask = pad_to_even_square(mp.P1.full.mask); % force to be square
        mp.P1.full.Narr = length(mp.P1.full.mask);
    end

    % Full model coordinates
    if(strcmpi(mp.centering,'interpixel') )
        mp.P2.full.xsDL = (- (mp.P1.full.Narr-1)/2:(mp.P1.full.Narr-1)/2)*mp.P2.full.dx/mp.P2.D;
    else
        mp.P2.full.xsDL = ( -mp.P1.full.Narr/2:(mp.P1.full.Narr/2 -1) )*mp.P2.full.dx/mp.P2.D;
    end
    [mp.P2.full.XsDL, mp.P2.full.YsDL] = meshgrid(mp.P2.full.xsDL);

end
