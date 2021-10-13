%% Test input checks
%
% Unit tests of the methods in the Check class.
%
classdef TestCheckTwoDimSquareArray < matlab.unittest.TestCase
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
    
%% Unit tests of two_dim_square_array

    methods (Test)
        
        function testNonNumeric(testCase)
            var = 'asdf';
            identifier = 'ValueError:InputMustBeNumeric';
            verifyError(testCase, @() Check.two_dim_square_array(var), identifier)
        end
        function testEmpty(testCase)
            var = [];
            identifier = 'ValueError:InputMustBeArray';
            verifyError(testCase, @() Check.two_dim_square_array(var), identifier)
        end
        function testScalar(testCase)
            var = 2.5;
            identifier = 'ValueError:InputMustBeArray';
            verifyError(testCase, @() Check.two_dim_square_array(var), identifier)
        end
        function testTwoDim(testCase)
            var = ones(3, 2, 4);
            identifier = 'ValueError:InputMustBeTwoDim';
            verifyError(testCase, @() Check.two_dim_square_array(var), identifier)
        end
        function testSquare(testCase)
            var = ones(3, 2);
            identifier = 'ValueError:InputMustBeSquare';
            verifyError(testCase, @() Check.two_dim_square_array(var), identifier)
        end
        
    end 

end
