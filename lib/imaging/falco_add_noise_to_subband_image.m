% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Add noise (photon shot, dark current, & read) to a simulated image.
%
% INPUTS
% ------
% mp : structure of model parameters
% imageIn : noiseless starting image for a given subband [normalized intensity]
% iSubband : index of subband in which the image was taken
%
% OUTPUTS
% -------
% imageOut : noisy image [normalized intensity]

function imageOut = falco_add_noise_to_subband_image(mp, imageIn, iSubband)

    peakCounts = mp.detector.peakFluxVec(iSubband) * mp.detector.tExpVec(iSubband);
    peakElectrons = mp.detector.gain * peakCounts;

    imageInElectrons = peakElectrons * imageIn;
    factorNeededInMatlab = 1e12;
    
    imageInCounts = 0;
    for iExp = 1:mp.detector.Nexp
    
        % Add photon shot noise
        noisyImageInElectrons = factorNeededInMatlab*imnoise(imageInElectrons/factorNeededInMatlab, 'poisson');

        % Compute dark current
        darkCurrent = mp.detector.darkCurrentRate * mp.detector.tExpVec(iSubband)*ones(size(imageIn));
        darkCurrent = factorNeededInMatlab*imnoise(darkCurrent/factorNeededInMatlab, 'poisson');

        % Compute Gaussian read noise
        readNoise = mp.detector.readNoiseStd * randn(size(imageIn));
        
        % Convert back from e- to counts and then discretize
        imageInCounts = imageInCounts + round((noisyImageInElectrons + darkCurrent + readNoise)/mp.detector.gain)/mp.detector.Nexp;
        
    end

    % Convert back from counts to normalized intensity
    imageOut = imageInCounts / peakCounts; 
            
end
