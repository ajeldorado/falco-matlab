%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
% Test falco_add_noise_to_subband_image.m
%
classdef TestNoisyImage < matlab.unittest.TestCase  
%% Properties
%
% A presaved file with FALCO parameters was saved and is loaded to be used
% by methods.
%     properties
%         mp=Parameters();
%     end

%% Setup and Teardown Methods
%
%  Add and remove path to library functions to be tested.

    methods (TestClassSetup)
        function addPath(testCase)
            pathToFalco = fileparts(fileparts(fileparts(mfilename('fullpath')))); % falco-matlab directory;
            addpath(genpath([pathToFalco filesep 'lib/imaging']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            pathToFalco = fileparts(fileparts(fileparts(mfilename('fullpath')))); % falco-matlab directory;
            rmpath(genpath([pathToFalco filesep 'lib/imaging']));
        end
    end

%% Tests

    methods (Test)   
 
        function testPhotonShotNoise(testCase)
            
            trueMean = 1e-4; % [normalized intensity]
            imageIn = trueMean * ones(1000, 1000);
            peakFluxCoef = 1e10; % [counts/pixel/second]
            readNoiseStd = 0; 
            tExp = 38;
            gain = 4.2;
            darkCurrentRate = 0; %[e-/pixel/second]
            shotNoiseStd = sqrt(gain*tExp*trueMean*peakFluxCoef);
            
            mp.Nsbp = 3;
            iSubband = 3;
            mp.detector.gain = gain; %1.0; % [e-/count]
            mp.detector.darkCurrentRate = darkCurrentRate; % [e-/pixel/second]
            mp.detector.readNoiseStd = readNoiseStd; %1.7; % [e-/count]
            mp.detector.peakFluxVec = peakFluxCoef * ones(mp.Nsbp, 1); % [counts/pixel/second]
            mp.detector.tExpVec = tExp * ones(mp.Nsbp, 1); % [seconds]

            imageOut = falco_add_noise_to_subband_image(mp, imageIn, iSubband);
            
            imageStd = std(imageOut(:) - trueMean)*mp.detector.gain*peakFluxCoef*tExp;
            fprintf('Expected photon shot noise = %.5e\n', shotNoiseStd);
            fprintf('Meas photon shot noise     = %.5e\n', imageStd);
            
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            testCase.verifyThat(imageStd, IsEqualTo(shotNoiseStd, 'Within', RelativeTolerance(1e-2)))
        end        
        
        function testThatMeanStaysSame(testCase)
            
            trueMean = 1.123e-3; % [normalized intensity]
            imageIn = trueMean * ones(1000, 1000);
            
            mp.Nsbp = 3;
            iSubband = 3;
            mp.detector.gain = 4; %1.0; % [e-/count]
            mp.detector.darkCurrentRate = 0.1;%0.015; % [e-/pixel/second]
            mp.detector.readNoiseStd = 5; %1.7; % [e-/count]
            mp.detector.peakFluxVec = 1e10 * ones(mp.Nsbp, 1); % [counts/pixel/second]
            mp.detector.tExpVec = 10.0 * ones(mp.Nsbp, 1); % [seconds]

            imageOut = falco_add_noise_to_subband_image(mp, imageIn, iSubband);
            
            imageMean = mean(imageOut(:));
            fprintf('Mean before noise = %.5e\n', trueMean);
            fprintf('Mean after noise  = %.5e\n', imageMean);
            
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            testCase.verifyThat(imageMean, IsEqualTo(trueMean, 'Within', RelativeTolerance(1e-3)))
        end
        
        function testReadNoiseStd(testCase)
            
            trueMean = 0;%1e-10; % [normalized intensity]
            imageIn = trueMean * ones(1000, 1000);
            peakFluxCoef = 0.3; % [counts/pixel/second]
            readNoiseStd = 55; 
            tExp = 0.1;
            
            mp.Nsbp = 3;
            iSubband = 3;
            mp.detector.gain = 4; %1.0; % [e-/count]
            mp.detector.darkCurrentRate = 0;%0.015; % [e-/pixel/second]
            mp.detector.readNoiseStd = readNoiseStd; %1.7; % [e-/count]
            mp.detector.peakFluxVec = peakFluxCoef * ones(mp.Nsbp, 1); % [counts/pixel/second]
            mp.detector.tExpVec = tExp * ones(mp.Nsbp, 1); % [seconds]

            imageOut = falco_add_noise_to_subband_image(mp, imageIn, iSubband);
            
            imageStd = std(imageOut(:) - trueMean)*mp.detector.gain*peakFluxCoef*tExp;
            fprintf('True read noise std = %.5e\n', readNoiseStd);
            fprintf('Meas read noise std = %.5e\n', imageStd);
            
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            testCase.verifyThat(imageStd, IsEqualTo(readNoiseStd, 'Within', RelativeTolerance(1e-2)))
        end
         
        function testDarkCurrent(testCase)
            
            trueMean = 0;%1e-15; % [normalized intensity]
            imageIn = trueMean * ones(1000, 1000);
            peakFluxCoef = 1e-3; % [counts/pixel/second]
            readNoiseStd = 0; 
            tExp = 1e4;
            darkCurrentRate = 5.3; %[e-/pixel/second]
            darkCurrentStd = sqrt(tExp*darkCurrentRate);
            
            mp.Nsbp = 3;
            iSubband = 3;
            mp.detector.gain = 4; %1.0; % [e-/count]
            mp.detector.darkCurrentRate = darkCurrentRate; % [e-/pixel/second]
            mp.detector.readNoiseStd = readNoiseStd; %1.7; % [e-/count]
            mp.detector.peakFluxVec = peakFluxCoef * ones(mp.Nsbp, 1); % [counts/pixel/second]
            mp.detector.tExpVec = tExp * ones(mp.Nsbp, 1); % [seconds]

            imageOut = falco_add_noise_to_subband_image(mp, imageIn, iSubband);
            
            imageStd = std(imageOut(:) - trueMean)*mp.detector.gain*peakFluxCoef*tExp;
            fprintf('Expected dark current noise = %.5e\n', darkCurrentStd);
            fprintf('Meas dark current noise     = %.5e\n', imageStd);
            
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            testCase.verifyThat(imageStd, IsEqualTo(darkCurrentStd, 'Within', RelativeTolerance(1e-2)))
        end
        
    end
    
end
