% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function OUT = propcustom_mft_PtoFtoP(IN, FPM, charge, apRad, inVal, outVal, useGPU ) 
%
% Propagate from the pupil plane before a FPM to the pupil plane after it.
%
% FPM must be an array with dimensions MxM, where M = lambdaOverD*(beam diameter)
%
% Optional arguments allow for an opaque spot to be placed on the mask to
% block the phase singularity. The spot can also be offset from the mask center.
%
% INPUTS
% ------
% IN: 2-D E-field at the pupil from which to propagate
% FPM: 2-D, complex-valued array representing the FPM to use.
% apRad: radius of the beam at the starting pupil. Units of pixels.
% inVal: radial separation in lambda/D at which to start the Tukey window
%        transition between coarse and fine resolution.
% outVal: radial separation in lambda/D at which to end the Tukey window
%         transition between coarse and fine resolution.
% useGPU: boolean flag whether to use a GPU to speed up calculations
% diamSpotLamD: (optional) diameter of the opaque spot in lambda0/D.
% offsetsLamD: (optional) the (x, y) offsets of the opaque spot from the
%              mask center in units of lambda0/D.
%
% OUTPUTS
% -------
% OUT: 2-D E-field at the next pupil plane after propagating through the
%      FPM.
%
% WARNINGS
% --------
% 1. Use this only with azimuthal-only masks. Radial variation will cause
%    problems because the same FPM array is used at both the coarse and
%    fine resolutions.
% 2. The sampling of the coarse DFT is set as the width of the FPM array
%    divided by the beam diameter in pixels.

function OUT = propcustom_mft_PtoFtoP(IN, FPM, apRad, inVal, outVal, useGPU, varargin)

    % Optional arguments
    diamSpotLamD = 0; % Spot diameter [lambda0/D]
    offsetsLamD = [0, 0]; % (x, y) offsets of spot [lambda0/D]
    if length(varargin) == 1
        diamSpotLamD = varargin{1};
    elseif length(varargin) == 2
        diamSpotLamD = varargin{1};
        offsetsLamD = varargin{2};
    elseif length(varargin) > 1
        error('Too many inputs.')
    end

    showPlots2debug = false; 

    D = 2*apRad;    
    [NA, ~] = size(IN);
    NB = size(FPM, 1);
    pixPerLamD = NB/D; % [samples per lambda/D at coarse resolution]
    
    [X,Y] = meshgrid(-NB/2:NB/2-1);
    RHO = sqrt(X.^2 + Y.^2);
   
    windowKnee = 1-inVal/outVal;
    
    windowMASK1 = falco_gen_Tukey4vortex(2*outVal*pixPerLamD, RHO, windowKnee) ;
    windowMASK2 = falco_gen_Tukey4vortex(NB, RHO, windowKnee) ;

    % DFT vectors 
    x = (-NA/2:NA/2-1)/D;
    u1 = (-NB/2:NB/2-1)/pixPerLamD;
    u2 = (-NB/2:NB/2-1)*2*outVal/NB;
        
    if useGPU
        IN = gpuArray(IN);
        x = gpuArray(x);
        u1 = gpuArray(u1);
        u2 = gpuArray(u2);
        windowMASK1 = gpuArray(windowMASK1);
        windowMASK2 = gpuArray(windowMASK2);
    end

    if showPlots2debug; figure; imagesc(abs(IN)); axis image; colorbar; title('pupil'); end

    %% Full FOV, lower resolution DFT
    
    % Generate low-resolution central opaque spot
    if diamSpotLamD > 0
        inputs.pixresFPM = pixres; %--pixels per lambda/D
        inputs.rhoInner = diamSpotLamD/2; % radius of inner FPM amplitude spot (in lambda_c/D)
        inputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
        inputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
        inputs.centering = 'pixel';
        inputs.xOffset = offsetsLamD(1);
        inputs.yOffset = offsetsLamD(2);
        spotCoarse = falco_gen_annular_FPM(inputs);
        spotCoarse = pad_crop(spotCoarse, NB, 'extrapval', 1);
    else
        spotCoarse = 1;
    end
    
    FP1 = 1/(1*D*pixPerLamD)*exp(-1i*2*pi*u1'*x)*IN*exp(-1i*2*pi*x'*u1); 
    if showPlots2debug; figure; imagesc(log10(abs(FP1).^2)); axis image; colorbar; title('Large scale DFT'); end

    LP1 = 1/(1*D*pixPerLamD)*exp(-1i*2*pi*x'*u1)*(FP1.*spotCoarse.*FPM.*(1-windowMASK1))*exp(-1i*2*pi*u1'*x);
    if showPlots2debug; figure; imagesc(abs(FP1.*(1-windowMASK1))); axis image; colorbar; title('Large scale DFT (windowed)'); end
    
    %% Smaller FOV, finer sampling DFT
    
    % Generate high-resolution central opaque spot
    if diamSpotLamD > 0
        inputs.pixresFPM = NB/(2*outVal); %--pixels per lambda/D
        inputs.rhoInner = diamSpotLamD/2; % radius of inner FPM amplitude spot (in lambda_c/D)
        inputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
        inputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
        inputs.centering = 'pixel';
        inputs.xOffset = offsetsLamD(1);
        inputs.yOffset = offsetsLamD(2);
        spotFine = falco_gen_annular_FPM(inputs);
        spotFine = pad_crop(spotFine, NB, 'extrapval', 1);
    else
        spotFine = 1;
    end
    FP2 = 2*outVal/(1*D*NB)*exp(-1i*2*pi*u2'*x)*IN*exp(-1i*2*pi*x'*u2); 
    if showPlots2debug; figure; imagesc(log10(abs(FP2).^2)); axis image; colorbar; title('Fine sampled DFT'); end

    LP2 = 2*outVal/(1*D*NB)*exp(-1i*2*pi*x'*u2)*(FP2.*spotFine.*FPM.*windowMASK2)*exp(-1i*2*pi*u2'*x);        
    if showPlots2debug; figure; imagesc(abs(FP2.*windowMASK2)); axis image; colorbar; title('Fine sampled DFT (windowed)'); end
    OUT = LP1 + LP2;
    if showPlots2debug; figure; imagesc(abs(OUT)); axis image; colorbar; title('Lyot plane'); end

    if useGPU
        OUT = gather(OUT);
    end
    
end
