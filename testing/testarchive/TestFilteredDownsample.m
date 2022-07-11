%% Test falco_filtered_downsample()
%
% Unit tests of falco_filtered_downsample().
%
classdef TestFilteredDownsample < matlab.unittest.TestCase
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
    
    %% Unit tests

    methods (Test)

        function testPixelCenteredArray(testCase)
            
            arrayIn = ones(3);
            arrayIn = pad_crop(arrayIn, 9);
            arrayIn = pad_crop(arrayIn, 15, 'extrapval', 1);
            arrayIn = pad_crop(arrayIn, 21);

            downSampFac = 1/3;
            centering = 'pixel';

            arrayOut = falco_filtered_downsample(arrayIn, downSampFac, centering);

            arrayExpected = pad_crop(1, 3);
            arrayExpected = pad_crop(arrayExpected, 5, 'extrapval', 1);
            arrayExpected = pad_crop(arrayExpected, 7);
            
            maxAbsDiff = max(max(abs(arrayExpected - arrayOut)));
            
            testCase.verifyTrue(maxAbsDiff < 10*eps)
        end
        
        function testInterpixelCenteredArraySymmetry(testCase)
            
            arrayIn = ones(8);
            arrayIn = pad_crop(arrayIn, 16);
            arrayIn = pad_crop(arrayIn, 24, 'extrapval', 1);
            arrayIn = pad_crop(arrayIn, 32);

            downSampFac = 1/4;
            centering = 'interpixel';

            arrayOut = falco_filtered_downsample(arrayIn, downSampFac, centering);
            
            maxAbsDiff = max(max(abs(arrayOut - rot90(arrayOut, 2))));
            
            testCase.verifyTrue(maxAbsDiff < 10*eps)
        end
                
    end 

end
