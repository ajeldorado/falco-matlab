%% Test input checks
%
% Unit tests of the functions in the ../../lib/check/ directory.
%
classdef TestCheckRealScalar < matlab.unittest.TestCase
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
    
%% Unit tests of check_real_scalar

    methods (Test)
        
        function testNonNumeric(testCase)
            var = 'asdf';
            identifier = 'check_real_scalar:InputMustBeNumeric';
            verifyError(testCase, @() check_real_scalar(var), identifier)
        end
        function testEmpty(testCase)
            var = [];
            identifier = 'check_real_scalar:InputMustBeScalar';
            verifyError(testCase, @() check_real_scalar(var), identifier)
        end
        function testArray(testCase)
            var = eye(2);
            identifier = 'check_real_scalar:InputMustBeScalar';
            verifyError(testCase, @() check_real_scalar(var), identifier)
        end
        function testFinite(testCase)
            var = -Inf;
            identifier = 'check_real_scalar:InputMustBeFinite';
            verifyError(testCase, @() check_real_scalar(var), identifier)
        end
        function testReal(testCase)
            var = 1 + 2j;
            identifier = 'check_real_scalar:InputMustBeReal';
            verifyError(testCase, @() check_real_scalar(var), identifier)
        end
        
    end 

end
