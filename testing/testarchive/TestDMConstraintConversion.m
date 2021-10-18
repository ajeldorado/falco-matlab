%% Test conversion between DM constraint specification formats
%
% 
%
classdef TestDMConstraintConversion < matlab.unittest.TestCase
%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    properties
        dummy = 1;
    end

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
                        
        function test_convert_matrix_into_pairs(testCase)
            
            tieMat = zeros(6, 6);
            tieMat(5:6, 5:6) = 1; % 4 tied actuators --> 3 tied pairs
            tieMat(2:3, 2) = 2;  % 2 tied actuators --> 1 tied pair
            tieMat(4, 4) = 3; % tied block with just one actuator; functionally does nothing different
            tieMat([5, 21, 27]) = 4; % 3 tied actuators --> 2 tied pairs
            tieMat(1, 2) = -1; % dead
%             figure(1); imagesc(tieMat); axis xy equal tight; colorbar;
%             figure(2); imagesc(reshape(1:36, [6, 6])); axis xy equal tight; colorbar;
            
            tiePairsOut = falco_convert_dm_tie_matrix_into_pairs(tieMat);
            
            tiePairsExpected = [...
                29, 30;
                29, 35;
                29, 36;
                8, 9;
                5, 21;
                5, 27;
                ];
            
            testCase.verifyTrue(isequal(tiePairsOut, tiePairsExpected))
        end
        

        function test_convert_pairs_into_matrix(testCase)
            
            tiePairs = [...
                29, 30;
                29, 35;
                7, 7; % should be ignored
                36, 35;
                8, 9;
                5, 21;
                5, 27;
                18, 36;
                34, 18;
                ];

            nact = 6;
            tieMatExpected = zeros(nact, nact);
            tieMatExpected(5:6, 5:6) = 3;
            tieMatExpected([18, 34]) = 3;
            tieMatExpected(2:3, 2) = 2;  % 2 tied actuators --> 1 tied pair
            tieMatExpected([5, 21, 27]) = 1; % 3 tied actuators --> 2 tied pairs
            
            tieMatOut = falco_convert_dm_tie_pairs_into_matrix(tiePairs, nact);
%             figure(3); imagesc(tieMatExpected); axis xy equal tight; colorbar;
%             figure(4); imagesc(tieMatOut); axis xy equal tight; colorbar;
            
            testCase.verifyTrue(isequal(tieMatOut, tieMatExpected))
        end
        
        
    end
end
