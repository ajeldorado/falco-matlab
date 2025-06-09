% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% 
% Get the multiplicative E-field to shift the final image via a tip-tilt at a pupil plane.
%
% INPUTS
% ------
% mp : structure of model parameters
% EP4 : complex, 2-D E-field at exit pupil plane before tip/tilt is applied.
% lambda : wavelength in meters
%
% OUTPUTS
% -------
% EP4 : complex, 2-D E-field at exit pupil plane after tip/tilt is applied.
%
% NOTES
% -----
% N/A

function EP4mult = falco_get_full_model_detector_offsets(mp, modvar)

    if length(mp.Fend.full.xiOffset) ~= length(mp.Fend.full.etaOffset)
        error("mp.Fend.full.xiOffset and mp.Fend.full.etaOffset must have same length.")
    end

    %--Set the wavelength
    if ~isempty(modvar.lambda) %--For FALCO or for evaluation without WFSC
        lambda = modvar.lambda;
        if length(mp.Fend.full.xiOffset) == 1
            xiOffset = mp.Fend.full.xiOffset;
            etaOffset = mp.Fend.full.etaOffset;
        else
            error('An arbitrary wavelength cannot be used when the downstream image shift is given by a vector corresponding to exact wavelength values.')
        end

    else

        lambda = mp.full.lambdasMat(modvar.sbpIndex, modvar.wpsbpIndex);

        if length(mp.Fend.full.xiOffset) == 1
            xiOffset = mp.Fend.full.xiOffset;
            etaOffset = mp.Fend.full.etaOffset;
        elseif length(mp.Fend.full.xiOffset) == mp.Nsbp
            iSubband = modvar.sbpIndex;
            xiOffset = mp.Fend.full.xiOffset(iSubband);
            etaOffset = mp.Fend.full.etaOffset(iSubband);
        elseif length(mp.Fend.full.xiOffset) == length(mp.full.lambdas)
            iLambda = find(mp.full.lambdas==lambda);
            xiOffset = mp.Fend.full.xiOffset(iLambda);
            etaOffset = mp.Fend.full.etaOffset(iLambda);        
        else
            error('mp.Fend.full.xiOffset must either have length 1 or mp.full.NlamUnique or mp.Nsbp.')
        end
    
    end

    % display(xiOffset)
    % display(etaOffset)
    if (xiOffset ~= 0) || (etaOffset ~= 0)
        TTphase = (2*pi*(xiOffset*mp.P4.full.XsDL + etaOffset*mp.P4.full.YsDL));
        EP4mult = exp(1j*TTphase*mp.lambda0/lambda);
    else
        EP4mult = [];
    end

end
