%% Test User Defined pad_crop.m Function
%
% 
classdef TestPadCrop < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
%     properties
%         mp=DefaultParameters();
%     end
    
%% Setup and Teardown Methods at Class Level
%
%  Add and remove path to utility functions directory where pad_crop.m resides.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath(['../../../falco-matlab/lib/utils']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath(['../../../falco-matlab/lib/utils']))
        end
    end

%% Tests
%
%  Creates tests:
%
% # *testScalarInput1:* verify scalar input (1,3). The constrint
%       verifyEqual verifies that the actual and expected values have the
%       same size and that the corresponding values are equal.
% # *testScalarInput2:* verify scalar input (1,[2 3]). The constrint
%       verifyEqual verifies that the actual and expected values have the
%       same size and that the corresponding values are equal.
% # *testScalarInput3:* verify scalar input (1, [2;3]). The constrint
%       verifyEqual verifies that the actual and expected values have the
%       same size and that the corresponding values are equal.
% # *testScalarInput4:* verify scalar input (2, 3, 'extrapval', 1). The constrint
%       verifyEqual verifies that the actual and expected values have the
%       same size and that the corresponding values are equal.
    methods (Test)
        function testScalarInput1(testCase)
            actual = pad_crop(1, 3);
            expected = zeros(3,3); expected(2,2) = 1;
            testCase.verifyEqual(actual,expected)
        end
        function testScalarInput2(testCase)            
            actual = pad_crop(1, [2,3]);
            expected = zeros(2,3); expected(2,2) = 1;
            testCase.verifyEqual(actual,expected)
        end        
        function testScalarInput3(testCase)            
            actual = pad_crop(1, [2;3]);
            expected = zeros(2,3); expected(2,2) = 1;
            testCase.verifyEqual(actual,expected)
        end  
        function testScalarInput4(testCase)            
            actual = pad_crop(2, 3, 'extrapval', 1);
            expected = ones(3,3); expected(2,2) = 2;
            testCase.verifyEqual(actual,expected)
        end  
%%% Test 2-D Inputs
        function test2DInput1(testCase)            
            actual = pad_crop([1,2], 4);
            expected = zeros(4,4); expected(3,2:3) = [1,2];
            testCase.verifyEqual(actual,expected)
        end  
        function test2DInput2(testCase) 
            actual = pad_crop(eye(3), 4);
            expected = eye(4); expected(1,1) = 0;
            testCase.verifyEqual(actual,expected)
        end
        function test2DInput3(testCase) 
            actual = pad_crop(eye(3), [5,5]);
            expected = eye(4); expected(1,1) = 0; expected(5,5) = 0;
            testCase.verifyEqual(actual,expected)
        end
        function test2DInput4(testCase) 
            inp = eye(3);
            actual = pad_crop(inp, [5,5]);
            actual = pad_crop(actual, 3);
            expected = inp;
            testCase.verifyEqual(actual,expected)
        end
%%% Test 3-D and 4-D Inputs
        function test3D4DInput1(testCase) 
            actual = pad_crop(ones(3,2,4), 4);
            expected = zeros(4,4); expected(2:4, 2:3) = ones(3,2); expected = repmat(expected, [1,1,4]);
            testCase.verifyEqual(actual,expected)
        end        
        function test3D4DInput2(testCase) 
            actual = pad_crop(ones(3,2,4), [4,4]);
            expected = zeros(4,4); expected(2:4, 2:3) = ones(3,2); expected = repmat(expected, [1,1,4]);
            testCase.verifyEqual(actual,expected)
        end 
        function test3D4DInput3(testCase) 
            actual = pad_crop(ones(3,2,4,2), 4);
            expected = zeros(4,4); expected(2:4, 2:3) = ones(3,2); expected = repmat(expected, [1,1,4,2]);
            testCase.verifyEqual(actual,expected)
        end
        function test3D4DInput4(testCase) 
            actual = pad_crop(ones(3,2,4,2), [4,4]);
            expected = zeros(4,4); expected(2:4, 2:3) = ones(3,2); expected = repmat(expected, [1,1,4,2]);
            testCase.verifyEqual(actual,expected)
        end   
        function test3D4DInput5(testCase) 
            inp = ones(3,2,4,2);
            actual = pad_crop(inp, [4,4]);
            actual = pad_crop(actual, [3,2]);
            expected = inp;
            testCase.verifyEqual(actual,expected)
        end
        function test3D4DInput6(testCase)
            actual = pad_crop(ones(3,2,5,2), [4,4], 'extrapval', 1);
            expected = ones(4,4,5,2);
            testCase.verifyEqual(actual,expected)
        end
        function test3D4DInput7(testCase)
            actual = pad_crop(ones(3,2,5,2), [4;4], 'extrapval', 1);
            expected = ones(4,4,5,2);
            testCase.verifyEqual(actual,expected)
        end
    end    
end











