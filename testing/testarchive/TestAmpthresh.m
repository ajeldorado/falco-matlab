% function test_ampthresh_with_noise()

classdef TestAmpthresh < matlab.unittest.TestCase
%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath('../../lib/utils'));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath('../../lib/utils'));
        end
    end
    
%% Tests

    methods (Test)
        
        function test_ampthresh_with_noise(testCase)
            % Test that ampthresh recovers the pupil from a noisy image.
            pupil0 = zeros(100, 100);
            pupil0(30:70, 30:60) = 1;
            pupil = pupil0 + 0.1*rand(size(pupil0, 1), size(pupil0, 2));
            boolMask = ampthresh(pupil);
            testCase.verifyGreaterThan(sum(sum(boolMask == pupil0))/sum(sum(pupil0)), 0.99);  
        end
        
        function test_uniform_input(testCase)
            % Test that exception is raised for uniform input.
            testCase.verifyError(@()ampthresh(ones(10)), 'ampthresh:InputCannotBeUniform')
        end

        function testNonNumericInput(testCase)
            % Test that exception is raised for non-numeric value of nBins.
            testCase.verifyError(@()ampthresh(eye(5), 'string'), 'ampthresh:InputMustBeNumeric')
        end
        
    end    
end