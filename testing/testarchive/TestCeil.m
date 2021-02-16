%% Test ceil_even.m and ceil_odd.m
%
% Test applies known simple tests to ceil_even.m and ceil_odd.m to test
% results as we assign odd and even inputs to each function. Then the tests
% verify that the actual solution is equal to the expected known
% solution.This function does not enforce number type so there is not good
% reason to supply a test to verify the actual solution is of interger
% type.
%
classdef TestCeil < matlab.unittest.TestCase
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
% # *test ceil_odd* with even input
% # *test ceil_odd* with odd input
% # *test ceil_even* with even input
% # *test ceil_even* withg odd input
    methods (Test)
        function testCeilOdd1(testCase)
            actSolution = ceil_odd(56);
            expSolution = 57;
            testCase.verifyEqual(actSolution,expSolution);  
        end
        function testCeilOdd2(testCase)
            actSolution = ceil_odd(57);
            expSolution = 57;
            testCase.verifyEqual(actSolution,expSolution);  
        end
        function testCeilEven1(testCase)
            actSolution = ceil_even(56);
            expSolution = 56;
            testCase.verifyEqual(actSolution,expSolution);  
        end
        function testCeilEven2(testCase)
            actSolution = ceil_even(55);
            expSolution = 56;
            testCase.verifyEqual(actSolution,expSolution);  
        end  
    end    
end