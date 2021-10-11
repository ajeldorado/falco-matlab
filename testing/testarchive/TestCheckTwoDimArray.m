%% Test input checks
%
% Unit tests of the functions in the ../../lib/check/ directory.
%
classdef TestCheckTwoDimArray < matlab.unittest.TestCase
%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath('../../lib/check'));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath('../../lib/check'));
        end
    end
    
%% Unit tests of check_two_dim_array

    methods (Test)
        
        function testNonNumeric(testCase)
            var = 'asdf';
            identifier = 'check_two_dim_array:InputMustBeNumeric';
            verifyError(testCase, @() check_two_dim_array(var), identifier)
        end
        function testEmpty(testCase)
            var = [];
            identifier = 'check_two_dim_array:InputMustBeArray';
            verifyError(testCase, @() check_two_dim_array(var), identifier)
        end
        function testScalar(testCase)
            var = 2.5;
            identifier = 'check_two_dim_array:InputMustBeArray';
            verifyError(testCase, @() check_two_dim_array(var), identifier)
        end
        function testTwoDim(testCase)
            var = ones(3, 2, 4);
            identifier = 'check_two_dim_array:InputMustBeTwoDim';
            verifyError(testCase, @() check_two_dim_array(var), identifier)
        end
        
    end 

end
