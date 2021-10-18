%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test DM Apodized V Simple Design
%
% We test pupil Luvoir A (translation and rotation), and test pupil Luvoir
% A in Lyot stop mode.
classdef TestPupilLuvoirA < matlab.unittest.TestCase   
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
% # *testLuvoirATranslation* verify the pupil Luvoir A translation
%                            meet given constraints.
% # *testLuvoirARotation*  verify the pupil Luvoir A rotation
%                            meet given constraints.
% # *testLuvoirATransLS* verify the pupil Luvoir A in Lyot stop mode translation
%                            meet given constraints.
% # *testLuvoirARotLS*  verify the pupil Luvoir A in Lyot stop mode rotation
%                            meet given constraints
    methods (Test)    
        function testLuvoirATranslation(testCase)
            inputs.Nbeam = 200;
            pupil = falco_gen_pupil_LUVOIR_A_final(inputs);
            
            % Translation test
            inputs.xShear = -11/100;
            inputs.yShear = 19/100;
            pupilOffset = falco_gen_pupil_LUVOIR_A_final(inputs);
            diff = pad_crop(pupil, size(pupilOffset)) - circshift(pupilOffset, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))),1e-8)
        end
        function testLuvoirARotation(testCase)
            inputs.Nbeam = 200;
            pupil = falco_gen_pupil_LUVOIR_A_final(inputs);
            inputs.xShear = -11/100;
            inputs.yShear = 19/100;
            pupilOffset = falco_gen_pupil_LUVOIR_A_final(inputs);
            
            % Test rotation (and translation)
            inputs.clock_deg = 90;
            pupilRotOffset = falco_gen_pupil_LUVOIR_A_final(inputs);
            % figure(5); imagesc(pupilRotOffset); axis xy equal tight; colorbar; drawnow;
            pupilRot = zeros(size(pupil));
            pupilRot(2:end, 2:end) = rot90(pupil(2:end, 2:end),-1);
            diff = pad_crop(pupilRot, size(pupilOffset)) - circshift(pupilRotOffset, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))/sum(pupilRotOffset(:))),1e-3)
        end
        function testLuvoirATransLS(testCase)
            inputs.Nbeam = 200;
            inputs.flagLyot = true;
            inputs.ID = 0.30;
            inputs.OD = 0.80;
            inputs.clock_deg = 0;
            pupil = falco_gen_pupil_LUVOIR_A_final(inputs);

            % Translation test
            inputs.xShear = -11/100;
            inputs.yShear = 19/100;
            pupilOffset = falco_gen_pupil_LUVOIR_A_final(inputs);
            diff = pad_crop(pupil, size(pupilOffset)) - circshift(pupilOffset, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))),1e-8)
        end
        function testLuvoirARotLS(testCase)
            inputs.Nbeam = 200;
            inputs.flagLyot = true;
            inputs.ID = 0.30;
            inputs.OD = 0.80;
            inputs.clock_deg = 0;
            pupil = falco_gen_pupil_LUVOIR_A_final(inputs);
            inputs.xShear = -11/100;
            inputs.yShear = 19/100;
            pupilOffset = falco_gen_pupil_LUVOIR_A_final(inputs);

            % Test rotation (and translation)
            inputs.clock_deg = 90;
            pupilRotOffset = falco_gen_pupil_LUVOIR_A_final(inputs);
            pupilRot = zeros(size(pupil));
            pupilRot(2:end, 2:end) = rot90(pupil(2:end, 2:end),-1);
            diff = pad_crop(pupilRot, size(pupilOffset)) - circshift(pupilRotOffset, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))/sum(pupilRotOffset(:))), 1e-4)
        end
    end    
end