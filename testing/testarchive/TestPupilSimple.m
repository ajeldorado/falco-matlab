%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test Pupil Simple
%
% We test Pupil simple (area, translation and rotation).
classdef TestPupilSimple < matlab.unittest.TestCase    
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
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
            rmpath(genpath([pathToFalco filesep 'lib_external']));
        end
    end

%% Tests
%
%  Creates four tests:
%
% # *testPupilSimple* Here we group several verifications together to
%                     verify area, translation, and rotation of pupil.
%
    methods (Test)    
        function testPupilSimple(testCase)
            inputs.Nbeam = 100;
            inputs.Npad = 180;
            inputs.OD = 1.0;
            pupil = falco_gen_pupil_Simple(inputs);
            
            % Area test for open circle
            areaExpected = pi/4*(inputs.OD^2)*(inputs.Nbeam^2);
            area = sum(pupil(:));
            testCase.verifyLessThan(abs(areaExpected-area)/areaExpected*100, 0.01) % test that area error is < 0.01%
            
            % Area test for annulus
            inputs.ID = 0.20;
            pupil = falco_gen_pupil_Simple(inputs);
            areaExpected = pi/4*(inputs.OD^2 - inputs.ID^2)*(inputs.Nbeam^2);
            area = sum(pupil(:));
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            testCase.verifyThat(area, IsEqualTo(areaExpected,'Within', RelativeTolerance(0.0001)))
            
            % Include struts
            inputs.wStrut = 0.02;
            inputs.angStrut = 10 + [90, 180, 270, 315];
            pupil = falco_gen_pupil_Simple(inputs);
            
            % Translation test
            inputs.xShear = -11/100;
            inputs.yShear = 19/100;
            pupilOffset = falco_gen_pupil_Simple(inputs);
            diff = pad_crop(pupil, size(pupilOffset)) - circshift(pupilOffset, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))), 1e-8) 
            
            % Test rotation (and translation)
            inputs.clocking = 90;
            pupilRotOffset = falco_gen_pupil_Simple(inputs);
            pupilRot = zeros(size(pupil));
            pupilRot(2:end, 2:end) = rot90(pupil(2:end, 2:end),-1);
            diff = pad_crop(pupilRot, size(pupilOffset)) - circshift(pupilRotOffset, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))/sum(pupilRotOffset(:))), 1e-4)
        end
    end    
end