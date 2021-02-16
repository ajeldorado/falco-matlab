%% Test User Defined Sinc Function
%
% We apply usual tests such as xx= 0, pi/6, pi/4/, pi/2, pi, 3pi/2, etc.
% The test assumes that the user sets path to FALCO functionality but there
% is an extra test which verifies that the sinc.m is actually a file in the
% ../lib/utils in the FALCO path.
classdef TestSinc < matlab.unittest.TestCase
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
            addpath(genpath([testCase.mp.path.falco filesep 'lib/utils']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco filesep 'lib/utils']))
        end
    end  
     
%% Tests
%
%  Creates four tests:
%
% # *testSincpath:* verify sinc.m is in the path of FALCO
% # *testSinczero:* test sinc.m with zero input
% # *testSincPiosix:* test sinc.mwith pi/6 input
% # *testSincPiofour:* test sinc.m with pi/4 input
% # *testSincPiotwo:* test sinc.m with pi/2 input
% # *testSincPi:* test sinc.m with pi   input
    methods (Test)
        function testSincpath(testCase)
            import matlab.unittest.constraints.IsFile;
            act = '../../../falco-matlab/lib/utils/sinc.m';
            testCase.verifyThat(act,IsFile)
        end
        function testSinczero(testCase)
            xx=0;
            actSolution = sinc(xx);
            expSolution = 1;
            testCase.verifyEqual(actSolution,expSolution);
            testCase.verifyInstanceOf(actSolution, 'double');
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
        function testSincSize(testCase)
            xx=[pi/6 pi/4 pi/2 pi];
            actSolution = sinc(xx);
            expSolution = sin(pi*xx)./(pi*xx);
            testCase.verifyEqual(actSolution,expSolution);  
            testCase.verifySize(actSolution,size(xx));
        end
    end    
end