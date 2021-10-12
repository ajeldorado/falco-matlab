%% Test input checks
%
% Unit tests of the methods in the Check class.
%
classdef TestCheckRealScalar < matlab.unittest.TestCase
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
    
%% Unit tests of Check.real_scalar

    methods (Test)
                
        function testNonNumeric(testCase)
            var = 'asdf';
            identifier = 'real_scalar:InputMustBeNumeric';
            verifyError(testCase, @() Check.real_scalar(var), identifier)
        end
        function testEmpty(testCase)
            var = [];
            identifier = 'real_scalar:InputMustBeScalar';
            verifyError(testCase, @() Check.real_scalar(var), identifier)
        end
        function testArray(testCase)
            var = eye(2);
            identifier = 'real_scalar:InputMustBeScalar';
            verifyError(testCase, @() Check.real_scalar(var), identifier)
        end
        function testFinite(testCase)
            var = -Inf;
            identifier = 'real_scalar:InputMustBeFinite';
            verifyError(testCase, @() Check.real_scalar(var), identifier)
        end
        function testReal(testCase)
            var = 1 + 2j;
            identifier = 'real_scalar:InputMustBeReal';
            verifyError(testCase, @() Check.real_scalar(var), identifier)
        end
        
    end 

end
