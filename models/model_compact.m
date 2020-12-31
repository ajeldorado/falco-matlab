% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_compact(mp,   modvar)
% Blind model used by the estimator and controller
% Does not include unknown aberrations/errors that are in the full model.
%
% INPUTS
% mp : structure of model parameters
% modvar : structure of model variables
%
% OUTPUTS
% Eout : 2-D electric field at final focus
%
% NOTES
% modvar structure fields required for model_compact:
% - sbpIndex
% - whichSource

function [Eout, varargout] = model_compact(mp, modvar, varargin)

modvar.wpsbpIndex = 0; %--Dummy index since not needed in compact model

% Set default values of input parameters
normFac = mp.Fend.compact.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
flagEval = false; % flag to use a different (usually higher) resolution at final focal plane for evaluation
flagNewNorm = false;
%--Enable different arguments values by using varargin
icav = 0; % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case{'getnorm'}
            normFac = 0; 
            flagNewNorm = true;
        case {'normoff','unnorm','nonorm'} % Set to 0 when finding the normalization factor
            normFac = 1;
        case {'eval'} % Set to 0 when finding the normalization factor
            flagEval = true;
        otherwise
            error('model_compact: Unknown keyword: %s\n', varargin{icav});
    end
end

%--Normalization factor for compact evaluation model
if ~flagNewNorm && flagEval
    normFac = mp.Fend.eval.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
end

%--Set the wavelength
if(isfield(modvar,'lambda'))
    lambda = modvar.lambda;
else
    lambda = mp.sbp_centers(modvar.sbpIndex);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input E-fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Include the star position and weight in the starting wavefront
iStar = modvar.starIndex;
xiOffset = mp.compact.star.xiOffsetVec(iStar);
etaOffset = mp.compact.star.etaOffsetVec(iStar);
starWeight = mp.compact.star.weights(iStar);
TTphase = (-1)*(2*pi*(xiOffset*mp.P2.compact.XsDL + etaOffset*mp.P2.compact.YsDL));
Ett = exp(1i*TTphase*mp.lambda0/lambda);
Ein = sqrt(starWeight) * Ett .* mp.P1.compact.E(:, :, modvar.sbpIndex); 

if strcmpi(modvar.whichSource,'offaxis') %--Use for throughput calculations 
    TTphase = (-1)*(2*pi*(modvar.x_offset*mp.P2.compact.XsDL + modvar.y_offset*mp.P2.compact.YsDL));
    Ett = exp(1i*TTphase*mp.lambda0/lambda);
    Ein = Ett .* Ein;
end
% if strcmpi(modvar.whichSource,'offaxis') %--Use for throughput calculations 
%     TTphase = (-1)*(2*pi*(modvar.x_offset*mp.P2.compact.XsDL + modvar.y_offset*mp.P2.compact.YsDL));
%     Ett = exp(1i*TTphase*mp.lambda0/lambda);
%     Ein = Ett.*mp.P1.compact.E(:,:,modvar.sbpIndex);
% else %--Backward compatible with code without tip/tilt offsets in the Jacobian
%     Ein = mp.P1.compact.E(:, :, modvar.sbpIndex);  
% end

% if(isfield(mp,'ttx'))
%     %--Scale by lambda/lambda0 because ttx and tty are in lambda0/D
%     x_offset = mp.ttx(modvar.ttIndex)*(mp.lambda0/lambda);
%     y_offset = mp.tty(modvar.ttIndex)*(mp.lambda0/lambda);
% 
%     TTphase = (-1)*(2*pi*(x_offset*mp.P2.compact.XsDL + y_offset*mp.P2.compact.YsDL));
%     Ett = exp(1i*TTphase*mp.lambda0/lambda);
%     Ein = Ett.*mp.P1.compact.E(:,:,modvar.sbpIndex);  
% 
% elseif strcmpi(modvar.whichSource,'offaxis') %--Use for throughput calculations 
%     TTphase = (-1)*(2*pi*(modvar.x_offset*mp.P2.compact.XsDL + modvar.y_offset*mp.P2.compact.YsDL));
%     Ett = exp(1i*TTphase*mp.lambda0/lambda);
%     Ein = Ett.*mp.P1.compact.E(:,:,modvar.sbpIndex); 
%     
% else %--Backward compatible with code without tip/tilt offsets in the Jacobian
%     Ein = mp.P1.compact.E(:, :, modvar.sbpIndex);  
% end

%--Shift the source off-axis to compute the intensity normalization value.
%  This replaces the previous way of taking the FPM out in the optical model.
if normFac == 0
    source_x_offset = mp.source_x_offset_norm; %--source offset in lambda0/D for normalization
    source_y_offset = mp.source_y_offset_norm; %--source offset in lambda0/D for normalization
    TTphase = (-1)*(2*pi*(source_x_offset*mp.P2.compact.XsDL + source_y_offset*mp.P2.compact.YsDL));
    Ett = exp(1i*TTphase*mp.lambda0/lambda);
    Ein = Ett.*mp.P1.compact.E(:,:,modvar.sbpIndex); 
end

%--Apply a Zernike (in amplitude) at input pupil if specified
if ~isfield(modvar, 'zernIndex')
    modvar.zernIndex = 1;
end
%--Only used for Zernike sensitivity control, which requires the perfect 
% E-field of the differential Zernike term.
if modvar.zernIndex ~=1
    indsZnoll = modvar.zernIndex; %--Just send in 1 Zernike mode
    zernMat = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    zernMat = padOrCropEven(zernMat,mp.P1.compact.Narr);
    Ein = Ein.*zernMat*(2*pi*1i/lambda)*mp.jac.Zcoef(mp.jac.zerns==modvar.zernIndex);
end

%--Define what the complex-valued FPM is if the coronagraph is some type of HLC.
switch lower(mp.layout)
    case{'fourier'}
        switch upper(mp.coro) 
            case{'EHLC'} %--DMs, optional apodizer, extended FPM with metal and dielectric modulation and outer stop, and LS. Uses 1-part direct MFTs to/from FPM
                mp.FPM.mask = falco_gen_EHLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'compact'); %--Complex transmission map of the FPM.
            case{'HLC'} %--DMs, optional apodizer, FPM with optional metal and dielectric modulation, and LS. Uses Babinet's principle about FPM.
                mp.FPM.mask = falco_gen_HLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'compact'); %--Complex transmission map of the FPM.
        end
        
    case{'wfirst_phaseb_simple','wfirst_phaseb_proper','fpm_scale','proper'} %--Use compact model as the full model, and the general FALCO model as the compact model, or %--Use the actual Phase B compact model as the compact model.
        switch upper(mp.coro)     
            case{'HLC'}
                mp.FPM.mask = mp.compact.FPMcube(:, :, modvar.sbpIndex);
        end
end

%--Select which optical layout's compact model to use and get the output E-field
switch lower(mp.layout)
    case{'fourier'}
        if(mp.flagFiber)
            [Eout, Efiber] = model_compact_general(mp, lambda, Ein, normFac, flagEval);
            varargout{1} = Efiber;
        else
            Eout = model_compact_general(mp, lambda, Ein, normFac, flagEval);
        end
        
    case{'proper'}
        switch upper(mp.coro)
            case{'HLC'}
                    Eout = model_compact_scale(mp, lambda, Ein, normFac, flagEval);
            otherwise
                    Eout = model_compact_general(mp, lambda, Ein, normFac, flagEval);
        end
        
    case{'wfirst_phaseb_simple','wfirst_phaseb_proper'} %--Use compact model as the full model, and the general FALCO model as the compact model, or %--Use the actual Phase B compact model as the compact model.
        switch upper(mp.coro)
            case{'SPLC'}
                Eout = model_compact_general(mp, lambda, Ein, normFac, flagEval);
            case{'HLC'}
                Eout = model_compact_scale(mp, lambda, Ein, normFac, flagEval);
        end
        
    case{'fpm_scale'}
        switch upper(mp.coro)
            case{'HLC'}
                Eout = model_compact_scale(mp, lambda, Ein, normFac, flagEval);
        end
end
    
end %--END OF FUNCTION