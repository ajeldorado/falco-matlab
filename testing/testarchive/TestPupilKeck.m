%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test falco_gen_pupil_Keck.m
%
% Test that falco_gen_pupil_Keck.m generates apertures as expected.

classdef TestPupilKeck < matlab.unittest.TestCase  
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
            addpath(genpath([pathToFalco filesep 'lib']));
            addpath(genpath([pathToFalco filesep 'lib_external']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            pathToFalco = fileparts(fileparts(fileparts(mfilename('fullpath')))); % falco-matlab directory;
            rmpath(genpath([pathToFalco filesep 'lib']));
            addpath(genpath([pathToFalco filesep 'lib_external']));
        end
    end

%% Tests

    methods (Test)   
 
        function testRotation(testCase)
            inputs.Nbeam = 200; % number of points across the pupil diameter

            % Optional Inputs
            inputs.centering = 'pixel';%'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
            inputs.clocking = 65; % CCW rotation [degrees]
            inputs.xShear = 0.2; % [pupil diameters]
            inputs.yShear = -0.3; % [pupil diameters]

            pupil65 = falco_gen_pupil_Keck(inputs);

            inputs.clocking = 125; % CCW rotation [degrees]
            pupil125 = falco_gen_pupil_Keck(inputs);

            sumAbsDiff = sum(sum(abs(pupil125-pupil65)));

            testCase.verifyLessThan(sumAbsDiff, 1e-7) 
        end
        
        function testTranslation(testCase)
            inputs.Nbeam = 200; % number of points across the pupil diameter

            inputs.centering = 'pixel'; % 'interpixel' or 'pixel'
            inputs.clocking = 10; % CCW rotation [degrees]
            pupilCentered = falco_gen_pupil_Keck(inputs);

            inputs.xShear = 0.5; % [pupil diameters]
            inputs.yShear = -0.2; % [pupil diameters]
            pupilShifted = falco_gen_pupil_Keck(inputs);
            pupilRecentered = circshift(pupilShifted, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);

            pupilCentered = pad_crop(pupilCentered, size(pupilShifted));

            sumAbsDiff = sum(sum(abs(pupilRecentered-pupilCentered)));
            
            testCase.verifyLessThan(sumAbsDiff, 1e-7) 
        end
        
    end
    
end
