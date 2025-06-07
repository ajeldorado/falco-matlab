% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% 
% Apply a lateral shear to the final image via a tip-tilt at a pupil plane.
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

function EP4 = falco_apply_detector_offsets(mp, EP4, lambda, model_type)

    if strcmpi(model_type, 'compact') || strcmpi(model_type, 'Jacobian')

        if length(mp.Fend.compact.xiOffset) ~= length(mp.Fend.compact.etaOffset)
            error("mp.Fend.compact.xiOffset and mp.Fend.compact.etaOffset must have same length.")
        end
    
        if length(mp.Fend.compact.xiOffset) == 1
            xiOffset = mp.Fend.compact.xiOffset;
            etaOffset = mp.Fend.compact.etaOffset;
        elseif length(mp.Fend.compact.xiOffset) == mp.Nsbp
            iSubband = find(mp.sbp_centers==lambda);
            xiOffset = mp.Fend.compact.xiOffset(iSubband);
            etaOffset = mp.Fend.compact.etaOffset(iSubband);
        else
            error('mp.Fend.compact.xiOffset must either have length 1 or mp.Nsbp.')
        end
        
        if (xiOffset ~= 0) || (etaOffset ~= 0)
            TTphase = (2*pi*(xiOffset*mp.P4.compact.XsDL + etaOffset*mp.P4.compact.YsDL));
            EP4 = EP4 .* exp(1j*TTphase*mp.lambda0/lambda);
        end

    elseif strcmpi(model_type, 'full')

        if length(mp.Fend.full.xiOffset) ~= length(mp.Fend.full.etaOffset)
            error("mp.Fend.full.xiOffset and mp.Fend.full.etaOffset must have same length.")
        end
    
        if length(mp.Fend.full.xiOffset) == 1
            xiOffset = mp.Fend.full.xiOffset;
            etaOffset = mp.Fend.full.etaOffset;
        elseif length(mp.Fend.full.xiOffset) == length(mp.full.lambdas)
            iLambda = find(mp.full.lambdas==lambda);
            xiOffset = mp.Fend.full.xiOffset(iLambda);
            etaOffset = mp.Fend.full.etaOffset(iLambda);
        elseif length(mp.Fend.full.xiOffset) == mp.Nsbp
            iSubband = find(mp.sbp_centers==lambda);
            xiOffset = mp.Fend.full.xiOffset(iSubband);
            etaOffset = mp.Fend.full.etaOffset(iSubband);
        else
            error('mp.Fend.full.xiOffset must either have length 1 or mp.full.NlamUnique or mp.Nsbp.')
        end
        
        if (xiOffset ~= 0) || (etaOffset ~= 0)
            TTphase = (2*pi*(xiOffset*mp.P4.full.XsDL + etaOffset*mp.P4.full.YsDL));
            EP4 = EP4 .* exp(1j*TTphase*mp.lambda0/lambda);
        end

    else
        error("model_type must be 'compact', 'Jacobian', or 'full'.")
    end

end
