%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test DM Apodized V Simple Design
%
% We test pupil Luvoir B (translation and rotation).
classdef TestPupilLuvoirB < matlab.unittest.TestCase  
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
    properties
        mp=Parameters();
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath([testCase.mp.path.falco 'lib/utils']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco 'lib/utils']))
        end
    end    
    
%% Tests
%
%  Creates four tests:
%
% # *testLuvoirBTranslation* verify the pupil Luvoir B translation
%                            meet given constraints.
% # *testLuvoirBTranslationandRotation*  verify the pupil Luvoir B translation 
%                                        and rotation meet given constraints.
%
    methods (Test)    
        function testLuvoirBTranslation(testCase)
            inputs.Nbeam = 300;
            pupil = falco_gen_pupil_LUVOIR_B(inputs);
            
            % Translation test
            inputs.xShear = -11/100;
            inputs.yShear = 19/100;
            pupilOffset = falco_gen_pupil_LUVOIR_B(inputs);
            diff = pad_crop(pupil, size(pupilOffset)) - circshift(pupilOffset, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))),1e-8)
        end
        function testLuvoirBTranslationandRotation(testCase)
            inputs.Nbeam = 300;
            pupil = falco_gen_pupil_LUVOIR_B(inputs);
            % Translation test
            inputs.xShear = -11/100;
            inputs.yShear = 19/100;
            pupilOffset = falco_gen_pupil_LUVOIR_B(inputs);

            % Test rotation (and translation)
            inputs.clock_deg = 90;
            pupilRotOffset = falco_gen_pupil_LUVOIR_B(inputs);
            pupilRot = zeros(size(pupil));
            pupilRot(2:end, 2:end) = rot90(pupil(2:end, 2:end),-1);
            diff = pad_crop(pupilRot, size(pupilOffset)) - circshift(pupilRotOffset, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))/sum(pupilRotOffset(:))),1e-3)
        end
    end    
end