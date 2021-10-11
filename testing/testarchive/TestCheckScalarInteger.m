%% Test input checks
%
% Unit tests of the functions in the ../../lib/check/ directory.
%
classdef TestCheckScalarInteger < matlab.unittest.TestCase
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
    
%% Unit tests of check_scalar_integer

    methods (Test)
        
        function testNonNumeric(testCase)
            var = 'asdf';
            identifier = 'check_scalar_integer:InputMustBeNumeric';
            verifyError(testCase, @() check_scalar_integer(var), identifier)
        end
        function testEmpty(testCase)
            var = [];
            identifier = 'check_scalar_integer:InputMustBeScalar';
            verifyError(testCase, @() check_scalar_integer(var), identifier)
        end
        function testArray(testCase)
            var = eye(2);
            identifier = 'check_scalar_integer:InputMustBeScalar';
            verifyError(testCase, @() check_scalar_integer(var), identifier)
        end
        function testFinite(testCase)
            var = -Inf;
            identifier = 'check_scalar_integer:InputMustBeFinite';
            verifyError(testCase, @() check_scalar_integer(var), identifier)
        end
        function testReal(testCase)
            var = 1 + 2j;
            identifier = 'check_scalar_integer:InputMustBeReal';
            verifyError(testCase, @() check_scalar_integer(var), identifier)
        end
        function testInteger(testCase)
            var = 1.2;
            identifier = 'check_scalar_integer:InputMustBeIntegral';
            verifyError(testCase, @() check_scalar_integer(var), identifier)
        end
        
    end 

end
