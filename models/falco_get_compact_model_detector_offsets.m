% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% 
% Get the multiplicative E-field to shift the final image via a tip-tilt at
% a pupil plane in the compact model.
%
% INPUTS
% ------
% mp : structure of model parameters
% modvar : an instance of the ModelVariables class
%
% OUTPUTS
% -------
% EP4mult : change in E-field at exit pupil plane after tip/tilt is
% applied. If no tip/tilt is applied, it is the scalar 1.
%
% NOTES
% -----
% N/A

function EP4mult = falco_get_compact_model_detector_offsets(mp, modvar)

    if length(mp.Fend.compact.xiOffset) ~= length(mp.Fend.compact.etaOffset)
        error("mp.Fend.compact.xiOffset and mp.Fend.compact.etaOffset must have same length.")
    end

    %--Set the wavelength
    if ~isempty(modvar.lambda) %--For FALCO or for evaluation without WFSC
        lambda = modvar.lambda;
        if length(mp.Fend.compact.xiOffset) == 1
            xiOffset = mp.Fend.compact.xiOffset;
            etaOffset = mp.Fend.compact.etaOffset;
        else
            error('An arbitrary wavelength cannot be used when the downstream image shift is given by a vector corresponding to exact wavelength values.')
        end

    else

        lambda = mp.sbp_centers(modvar.sbpIndex);

        if length(mp.Fend.compact.xiOffset) == 1
            xiOffset = mp.Fend.compact.xiOffset;
            etaOffset = mp.Fend.compact.etaOffset;
        elseif length(mp.Fend.compact.xiOffset) == mp.Nsbp
            iSubband = modvar.sbpIndex;
            xiOffset = mp.Fend.compact.xiOffset(iSubband);
            etaOffset = mp.Fend.compact.etaOffset(iSubband);      
        else
            error('mp.Fend.compact.xiOffset must either have length 1 or mp.Nsbp.')
        end
    
    end

    % 180-degree rotations
    % This is applied to account for the 180-degree rotation that may be
    % applied in the compact model.
    rot180fac = 1; % Default is no extra rotations
    if mp.flagRotation
        if mod(mp.NrelayFend, 2) == 1
            rot180fac = -1;
        end
    end

    if (xiOffset ~= 0) || (etaOffset ~= 0)
        TTphase = rot180fac*(2*pi*(xiOffset*mp.P4.compact.XsDL + etaOffset*mp.P4.compact.YsDL));
        EP4mult = exp(1j*TTphase*mp.lambda0/lambda);
    else
        EP4mult = 1;
    end

end
