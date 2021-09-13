% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Blind model used by the estimator and controller
% Does not include unknown aberrations/errors that are in the full model.
%
% INPUTS
% ------
% mp : structure of model parameters
% modvar : structure of model variables
%
% OUTPUTS
% -------
% Eout : 2-D electric field at final focus

function [Eout, varargout] = model_compact(mp, modvar, varargin)

modvar.wpsbpIndex = -1; %--Dummy index since not needed in compact model

% Set default values of input parameters
normFac = mp.Fend.compact.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
flagEval = false; % flag to use a different (usually higher) resolution at final focal plane for evaluation
flagNewNorm = false;
flagUseFPM = true; % default is to have the FPM in the beam
%--Enable different arguments values by using varargin
icav = 0; % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case{'getnorm'}
            normFac = 0; 
            flagNewNorm = true;
        case {'normoff', 'unnorm', 'nonorm'} 
            normFac = 1;
        case {'eval'}
            flagEval = true;
        case{'nofpm', 'unocculted'}
            flagUseFPM = false;
        otherwise
            error('model_compact: Unknown keyword: %s\n', varargin{icav});
    end
end

%--Normalization factor for compact evaluation model
if ~flagNewNorm && flagEval
    normFac = mp.Fend.eval.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
end

%--Set the wavelength
if isfield(modvar, 'lambda')
    lambda = modvar.lambda;
else
    lambda = mp.sbp_centers(modvar.sbpIndex);
end

%--Include the star position and weight in the starting wavefront
iStar = modvar.starIndex;
xiOffset = mp.compact.star.xiOffsetVec(iStar);
etaOffset = mp.compact.star.etaOffsetVec(iStar);
starWeight = mp.compact.star.weights(iStar);
TTphase = (-1)*(2*pi*(xiOffset*mp.P2.compact.XsDL + etaOffset*mp.P2.compact.YsDL));
Ett = exp(1j*TTphase*mp.lambda0/lambda);
Ein = sqrt(starWeight) * Ett .* mp.P1.compact.E(:, :, modvar.sbpIndex); 

if strcmpi(modvar.whichSource, 'offaxis') %--Use for throughput calculations 
    TTphase = (-1)*(2*pi*(modvar.x_offset*mp.P2.compact.XsDL + modvar.y_offset*mp.P2.compact.YsDL));
    Ett = exp(1j*TTphase*mp.lambda0/lambda);
    Ein = Ett .* Ein;
end

%--Shift the source off-axis to compute the intensity normalization value.
%  This replaces the previous way of taking the FPM out in the optical model.
if normFac == 0
    TTphase = (-1)*(2*pi*(mp.source_x_offset_norm*mp.P2.compact.XsDL + mp.source_y_offset_norm*mp.P2.compact.YsDL));
    Ett = exp(1j*TTphase*mp.lambda0/lambda);
    Ein = Ett .* mp.P1.compact.E(:, :, modvar.sbpIndex); 
end

%--Apply a Zernike (in amplitude) at input pupil if specified
if ~isfield(modvar, 'zernIndex')
    modvar.zernIndex = 1;
end
%--Only used for Zernike sensitivity control, which requires the perfect 
% E-field of the differential Zernike term.
if modvar.zernIndex ~= 1
    indsZnoll = modvar.zernIndex; %--Just send in 1 Zernike mode
    zernMat = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam, mp.centering, indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    zernMat = padOrCropEven(zernMat, mp.P1.compact.Narr);
    Ein = Ein .* zernMat * (2*pi*1j/lambda) * mp.jac.Zcoef(mp.jac.zerns == modvar.zernIndex);
end

%--Define what the complex-valued FPM is if the coronagraph is some type of HLC.
switch lower(mp.layout)
    case{'fourier'}
        switch upper(mp.coro) 
            case{'EHLC'} %--DMs, optional apodizer, extended FPM with metal and dielectric modulation and outer stop, and LS. Uses 1-part direct MFTs to/from FPM
                mp.FPM.mask = falco_gen_EHLC_FPM_complex_trans_mat(mp, modvar.sbpIndex, modvar.wpsbpIndex, 'compact'); %--Complex transmission map of the FPM.
            case{'HLC'} %--DMs, optional apodizer, FPM with optional metal and dielectric modulation, and LS. Uses Babinet's principle about FPM.
                mp.FPM.mask = falco_gen_HLC_FPM_complex_trans_mat(mp, modvar.sbpIndex, modvar.wpsbpIndex, 'compact'); %--Complex transmission map of the FPM.
        end
        
    case{'roman_phasec_proper', 'wfirst_phaseb_proper', 'fpm_scale', 'proper'}
        if strcmpi(mp.coro, 'HLC')
            mp.FPM.mask = mp.compact.FPMcube(:, :, modvar.sbpIndex);
        end
end

%--Select which optical layout's compact model to use and get the output E-field
if ~mp.flagFiber
    Eout = model_compact_general(mp, lambda, Ein, normFac, flagEval, flagUseFPM);
else
    [Eout, Efiber] = model_compact_general(mp, lambda, Ein, normFac, flagEval, flagUseFPM);
    varargout{1} = Efiber;
end
    
end %--END OF FUNCTION
