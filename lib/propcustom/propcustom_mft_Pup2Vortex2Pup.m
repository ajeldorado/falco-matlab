% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% Propagate from the pupil before a vortex to the pupil plane after it.
% An opaque spot is applied at the center of the vortex.
% 
% function OUT = propcustom_mft_Pup2Vortex2Pup( IN, charge, apRad,  inVal, outVal, useGPU, diamSpotLamD) 

function OUT = propcustom_mft_Pup2Vortex2Pup(IN, charge, apRad, inVal, outVal, useGPU, varargin)

    % Diameter of the central opaque spot is an optional input argument
    if length(varargin) == 1
        diamSpotLamD = varargin{1};
    elseif length(varargin) == 2
        diamSpotLamD = varargin{1};
        offsetsLamD = varargin{2};
    elseif length(varargin) > 1
        error('Too many inputs.')
    else
        diamSpotLamD = 0;
    end
    
    showPlots2debug = false; 

    D = 2*apRad;
    pixres = 4;%samples per lambda/D
    
    [NA, ~] = size(IN);
    NB = pixres*D; 
    
    [X, Y] = meshgrid(-NB/2:NB/2-1);
    RHO = sqrt(X.^2 + Y.^2);
   
    windowKnee = 1-inVal/outVal;
    
    windowMASK1 = falco_gen_Tukey4vortex( 2*outVal*pixres, RHO, windowKnee);
    windowMASK2 = falco_gen_Tukey4vortex( NB, RHO, windowKnee);

    % DFT vectors 
    x = (-NA/2:NA/2-1)/D;
    u1 = (-NB/2:NB/2-1)/pixres;
    u2 = (-NB/2:NB/2-1)*2*outVal/NB;
    
    FPM = falco_gen_vortex_mask(charge, NB);
    
    if(useGPU)
        IN = gpuArray(IN);
        x = gpuArray(x);
        u1 = gpuArray(u1);
        u2 = gpuArray(u2);
        windowMASK1 = gpuArray(windowMASK1);
        windowMASK2 = gpuArray(windowMASK2);
    end

    if showPlots2debug; figure;imagesc(abs(IN));axis image;colorbar; title('pupil'); end

    %% Large scale DFT

    % Generate central opaque spot
    if diamSpotLamD > 0
        inputs.pixresFPM = pixres; %--pixels per lambda/D
        inputs.rhoInner = diamSpotLamD/2; % radius of inner FPM amplitude spot (in lambda_c/D)
        inputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
        inputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
        inputs.centering = 'pixel';
        inputs.xOffset = offsetsLamD(1);
        inputs.yOffset = offsetsLamD(2);
        spot1 = falco_gen_annular_FPM(inputs);
        spot1 = pad_crop(spot1, NB, 'extrapval', 1);
    else
        spot1 = 1;
    end
    
    FP1 = 1/(1*D*pixres)*exp(-1i*2*pi*u1'*x)*IN*exp(-1i*2*pi*x'*u1); 
    if showPlots2debug; figure;imagesc(log10(abs(FP1).^2));axis image;colorbar; title('Large scale DFT'); end

    LP1 = 1/(1*D*pixres)*exp(-1i*2*pi*x'*u1)*(FP1.*spot1.*FPM.*(1-windowMASK1))*exp(-1i*2*pi*u1'*x);
    if showPlots2debug; figure;imagesc(abs(FP1.*(1-windowMASK1)));axis image;colorbar; title('Large scale DFT (windowed)'); end
    
    %% Fine sampled DFT

    % Generate central opaque spot
    if diamSpotLamD > 0
        inputs.pixresFPM = NB/(2*outVal); %--pixels per lambda/D
        inputs.rhoInner = diamSpotLamD/2; % radius of inner FPM amplitude spot (in lambda_c/D)
        inputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
        inputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
        inputs.centering = 'pixel';
        inputs.xOffset = offsetsLamD(1);
        inputs.yOffset = offsetsLamD(2);
        spot2 = falco_gen_annular_FPM(inputs);
        spot2 = pad_crop(spot2, NB, 'extrapval', 1);
    else
        spot2 = 1;
    end
    
    FP2 = 2*outVal/(1*D*NB)*exp(-1i*2*pi*u2'*x)*IN*exp(-1i*2*pi*x'*u2); 
    if showPlots2debug; figure;imagesc(log10(abs(FP2).^2));axis image;colorbar; title('Fine sampled DFT'); end
    FPM = falco_gen_vortex_mask(charge, NB);
    LP2 = 2*outVal/(1*D*NB)*exp(-1i*2*pi*x'*u2)*(FP2.*spot2.*FPM.*windowMASK2)*exp(-1i*2*pi*u2'*x);        
    if showPlots2debug; figure;imagesc(abs(FP2.*windowMASK2));axis image;colorbar; title('Fine sampled DFT (windowed)'); end
    OUT = LP1 + LP2;
    if showPlots2debug; figure;imagesc(abs(OUT));axis image;colorbar; title('Lyot plane'); end

    if(useGPU)
        OUT = gather(OUT);
    end
end