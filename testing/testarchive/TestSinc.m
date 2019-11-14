%% Test User Defined Sinc Function
%
% We apply usual tests such as xx= 0, pi/6, pi/4/, pi/2, pi, 3pi/2, etc.
% The test assumes that the function is not in the path and adds
% mp.path.falco to path for testing purposes, and then it removes path
% after tests are over. This is done for total independence.
classdef TestSinc < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
    properties
        param=load('Parameters.mat');
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods(TestMethodSetup)
        function addPath(testCase)
            addpath(genpath([testCase.param.mp.path.falco 'lib/utils']));
        end
    end
    methods(TestMethodTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.param.mp.path.falco 'lib/utils']))
        end
    end
    
%% Tests
%
%  Creates four testrs:
%
% # *test sinc* with zero input
% # *test sinc* with pi/6 input
% # *test sinc* with pi/4 input
% # *test sinc* with pi/2 input
% # *test sinc* with pi   input
    methods (Test)
        function testSinczero(testCase)
            xx=0;
            actSolution = sinc(xx);
            expSolution = 1;
            testCase.verifyEqual(actSolution,expSolution);  
        end
        function testSincPiosix(testCase)
            xx=pi/6;
            actSolution = sinc(xx);
            expSolution = sin(pi*xx)./(pi*xx);
            testCase.verifyEqual(actSolution,expSolution);  
        end
        function testSincPiofour(testCase)
            xx=pi/4;
            actSolution = sinc(xx);
            expSolution = sin(pi*xx)./(pi*xx);
            testCase.verifyEqual(actSolution,expSolution);  
        end
        function testSincPiotwo(testCase)
            xx=pi/2;
            actSolution = sinc(xx);
            expSolution = sin(pi*xx)./(pi*xx);
            testCase.verifyEqual(actSolution,expSolution);  
        end
        function testSincPi(testCase)
            xx=pi;
            actSolution = sinc(xx);
            expSolution = sin(pi*xx)./(pi*xx);
            testCase.verifyEqual(actSolution,expSolution);  
        end
    end    
end