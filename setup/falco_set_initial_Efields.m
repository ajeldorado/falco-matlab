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

    %% Entrance Pupil E-field

    % Full Model
    if ~isfield(mp.P1.full, 'E')

        mp.P1.full.E  = ones(mp.P1.full.Narr, mp.P1.full.Narr, mp.Nwpsbp, mp.Nsbp);

    else % if loading, pad to the correct size if necessary

        % if length(size(mp.P1.full.E)) ~= 4
        %     error('mp.P1.full.E must have 4 dimensions for y, x, wavelength, and subband.')
        % end

        if(size(mp.P1.full.E, 1) ~= mp.P1.full.Narr)
            EarrayTemp = mp.P1.full.E;
            mp.P1.full.E  = ones(mp.P1.full.Narr ,mp.P1.full.Narr, mp.Nwpsbp, mp.Nsbp); % make correctly sized array
            for iSubband = 1:mp.Nsbp
                for iLam = 1:mp.Nwpsbp
                    mp.P1.full.E(:, :, iLam, iSubband) = pad_crop(EarrayTemp(:, :, iLam, iSubband), mp.P1.full.Narr);
                end
            end
            clear EarrayTemp
        end

    end
    
    % Compact Model
    if ~isfield(mp.P1.compact, 'E')

        mp.P1.compact.E = ones(mp.P1.compact.Narr, mp.P1.compact.Narr, mp.Nsbp);

    else % if loading, pad to the correct size if necessary

        % if length(size(mp.P1.compact.E)) ~= 3
        %     error('mp.P1.compact.E must have 3 dimensions for y, x, and subband.')
        % end

        if(size(mp.P1.compact.E, 1) ~= mp.P1.compact.Narr)
            EcubeTemp = mp.P1.compact.E;
            mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp);
            for iSubband = 1:mp.Nsbp
                mp.P1.compact.E(:, :, iSubband) = pad_crop(EcubeTemp(:, :, iSubband), mp.P1.compact.Narr);
            end
            clear EcubeTemp
        end
        
    end

    %% Exit Pupil E-field Multiplier
    % If the array is already populated, it should not include tip/tilt.
    % Tip/tilt is added separately here.

    % Full Model
    if ~isfield(mp.P4.full, 'E')
        mp.P4.full.E = ones(mp.P4.full.Narr, mp.P4.full.Narr, mp.Nwpsbp, mp.Nsbp);
    else
        % if length(size(mp.P4.full.E)) ~= 4
        %     error('mp.P4.full.E must have 4 dimensions for y, x, wavelength, and subband.')
        % end
        
        if(size(mp.P4.full.E, 1) ~= mp.P4.full.Narr)
            EarrayTemp = mp.P4.full.E;
            mp.P4.full.E  = ones(mp.P4.full.Narr ,mp.P4.full.Narr, mp.Nwpsbp, mp.Nsbp);
            for iSubband = 1:mp.Nsbp
                for iLam = 1:mp.Nwpsbp
                    mp.P4.full.E(:, :, iLam, iSubband) = pad_crop(EarrayTemp(:, :, iLam, iSubband), mp.P4.full.Narr);
                end
            end
            clear EarrayTemp
        end

    end
    
    % Apply downstream tip/tilt
    modvar = ModelVariables();
    for iSubband = 1:mp.Nsbp
        modvar.sbpIndex = iSubband;
        for iLam = 1:mp.Nwpsbp
            modvar.wpsbpIndex = iLam;
            mp.P4.full.E(:, :, iLam, iSubband) = falco_get_full_model_detector_offsets(mp, modvar) .* mp.P4.full.E(:, :, iLam, iSubband);
        end
    end


    % Compact Model
    if ~isfield(mp.P4.compact, 'E')

        mp.P4.compact.E = ones(mp.P4.compact.Narr, mp.P4.compact.Narr, mp.Nsbp);

    else % if loading, pad to the correct size if necessary

        % if length(size(mp.P4.compact.E)) ~= 3
        %     error('mp.P4.compact.E must have 3 dimensions for y, x, and subband.')
        % end

        if(size(mp.P4.compact.E, 1) ~= mp.P4.compact.Narr)
            EcubeTemp = mp.P4.compact.E;
            mp.P4.compact.E = ones(mp.P4.compact.Narr,mp.P4.compact.Narr,mp.Nsbp);
            for iSubband = 1:mp.Nsbp
                mp.P4.compact.E(:, :, iSubband) = pad_crop(EcubeTemp(:, :, iSubband), mp.P4.compact.Narr);
            end
            clear EcubeTemp
        end
        
    end

    % Apply downstream tip/tilt
    modvar = ModelVariables();
    for iSubband = 1:mp.Nsbp
        modvar.sbpIndex = iSubband;
        mp.P4.compact.E(:, :, iSubband) = falco_get_compact_model_detector_offsets(mp, modvar) .* mp.P4.compact.E(:, :, iSubband);
    end

end %--END OF FUNCTION