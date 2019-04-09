% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function OUT = propcustom_mft_Pup2Vortex2Pup( IN, charge, apRad,  inVal, outVal, useGPU ) 
%
% Function to propagate from the pupil plane before a vortex FPM to the
% pupil plane after it.
%
% Written by Garreth Ruane.
% - Modified on 2019-04-05 by A.J. Riggs to remove the 1/1i term from each FT.

function OUT = propcustom_mft_Pup2Vortex2Pup( IN, charge, apRad,  inVal, outVal, useGPU )
%propcustom_mft_Pup2Vortex2Pup Propagates from the input pupil to output pupil

    showPlots2debug = false; 

    D = 2*apRad;
    lambdaOverD = 4;%samples per lambda/D
    
    [NA,~] = size(IN);
    NB = lambdaOverD*D; 
    
    [X,Y] = meshgrid(-NB/2:NB/2-1);
    RHO = sqrt(X.^2 + Y.^2);
   
    windowKnee = 1-inVal/outVal;
    
    windowMASK1 = falco_gen_Tukey4vortex( 2*outVal*lambdaOverD, RHO, windowKnee ) ;
    windowMASK2 = falco_gen_Tukey4vortex( NB, RHO, windowKnee ) ;

    % DFT vectors 
    x = (-NA/2:NA/2-1)/D;
    u1 = (-NB/2:NB/2-1)/lambdaOverD;
    u2 = (-NB/2:NB/2-1)*2*outVal/NB;
    
    FPM = falco_gen_vortex_mask( charge, NB );
    
    if(useGPU)
        IN = gpuArray(IN);
        x = gpuArray(x);
        u1 = gpuArray(u1);
        u2 = gpuArray(u2);
        windowMASK1 = gpuArray(windowMASK1);
        windowMASK2 = gpuArray(windowMASK2);
    end

    if showPlots2debug; figure;imagesc(abs(IN));axis image;colorbar; title('pupil'); end;

    %%%%%%% Large scale DFT

    FP1 = 1/(1*D*lambdaOverD)*exp(-1i*2*pi*u1'*x)*IN*exp(-1i*2*pi*x'*u1); 
    if showPlots2debug; figure;imagesc(log10(abs(FP1).^2));axis image;colorbar; title('Large scale DFT'); end;

    LP1 = 1/(1*D*lambdaOverD)*exp(-1i*2*pi*x'*u1)*(FP1.*FPM.*(1-windowMASK1))*exp(-1i*2*pi*u1'*x);
    if showPlots2debug; figure;imagesc(abs(FP1.*(1-windowMASK1)));axis image;colorbar; title('Large scale DFT (windowed)'); end;
    %%%%%%% Fine sampled DFT

    FP2 = 2*outVal/(1*D*NB)*exp(-1i*2*pi*u2'*x)*IN*exp(-1i*2*pi*x'*u2); 
    if showPlots2debug; figure;imagesc(log10(abs(FP2).^2));axis image;colorbar; title('Fine sampled DFT'); end;
    FPM = falco_gen_vortex_mask( charge, NB );
    LP2 = 2*outVal/(1*D*NB)*exp(-1i*2*pi*x'*u2)*(FP2.*FPM.*windowMASK2)*exp(-1i*2*pi*u2'*x);        
    if showPlots2debug; figure;imagesc(abs(FP2.*windowMASK2));axis image;colorbar; title('Fine sampled DFT (windowed)'); end;
    OUT = LP1 + LP2;
    if showPlots2debug; figure;imagesc(abs(OUT));axis image;colorbar; title('Lyot plane'); end;

    if(useGPU)
        OUT = gather(OUT);
    end
end