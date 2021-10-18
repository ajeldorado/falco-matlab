% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Define the initial stellar electric field(s) in the compact and full
% models.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function [mp] = falco_set_initial_Efields(mp)

    if ~isfield(mp.P1.full, 'E')
        mp.P1.full.E  = ones(mp.P1.full.Narr ,mp.P1.full.Narr, mp.Nwpsbp, mp.Nsbp);
    else % if loading, pad to the correct size if necessary
        if(size(mp.P1.full.E, 1) ~= mp.P1.full.Narr)
            EarrayTemp = mp.P1.full.E;
            mp.P1.full.E  = ones(mp.P1.full.Narr ,mp.P1.full.Narr, mp.Nwpsbp, mp.Nsbp);
            for si = 1:mp.Nsbp
                for wi = 1:mp.Nwpsbp
                    mp.P1.full.E(:, :, wi, si) = pad_crop(EarrayTemp(:, :, wi, si), mp.P1.full.Narr);
                end
            end
            clear EarrayTemp
        end
    end
    
    if ~isfield(mp.P1.compact, 'E')
        mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp);
    else % if loading, pad to the correct size if necessary
        if(size(mp.P1.compact.E, 1) ~= mp.P1.compact.Narr)
            EcubeTemp = mp.P1.compact.E;
            mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp);
            for si = 1:mp.Nsbp
                mp.P1.compact.E(:, :, si) = pad_crop(EcubeTemp(:, :, si), mp.P1.compact.Narr);
            end
            clear EcubeTemp
        end
    end

end %--END OF FUNCTION