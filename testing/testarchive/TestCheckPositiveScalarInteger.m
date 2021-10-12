%% Test input checks
%
% Unit tests of the methods in the Check class.
%
classdef TestCheckPositiveScalarInteger < matlab.unittest.TestCase
%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath('../../lib'));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath('../../lib'));
        end
    end
    
%% Unit tests of postive_positive_scalar_integer

    methods (Test)
        
        function testNonNumeric(testCase)
            var = 'asdf';
            identifier = 'positive_scalar_integer:InputMustBeNumeric';
            verifyError(testCase, @() Check.positive_scalar_integer(var), identifier)
        end
        function testEmpty(testCase)
            var = [];
            identifier = 'positive_scalar_integer:InputMustBeScalar';
            verifyError(testCase, @() Check.positive_scalar_integer(var), identifier)
        end
        function testArray(testCase)
            var = eye(2);
            identifier = 'positive_scalar_integer:InputMustBeScalar';
            verifyError(testCase, @() Check.positive_scalar_integer(var), identifier)
        end
        function testFinite(testCase)
            var = -Inf;
            identifier = 'positive_scalar_integer:InputMustBeFinite';
            verifyError(testCase, @() Check.positive_scalar_integer(var), identifier)
        end
        function testReal(testCase)
            var = 1 + 2j;
            identifier = 'positive_scalar_integer:InputMustBeReal';
            verifyError(testCase, @() Check.positive_scalar_integer(var), identifier)
        end
        function testInteger(testCase)
            var = 1.2;
            identifier = 'positive_scalar_integer:InputMustBeIntegral';
            verifyError(testCase, @() Check.positive_scalar_integer(var), identifier)
        end
        function testPositive(testCase)
            var = 0;
            identifier = 'positive_scalar_integer:InputMustBePositive';
            verifyError(testCase, @() Check.positive_scalar_integer(var), identifier)
        end
        
    end 

end
