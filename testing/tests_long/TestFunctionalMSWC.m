%% Functional Test MSWC
%
% The test script will perform the Wavefront Sensing and Control first, then it will
% test verify the values of the outputs to expected values define for each
% test.
classdef TestFunctionalMSWC < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we use the mp.path.falco to addpath to the
% function being tested.
    properties
        mp=ConfigurationMSWC();
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath([testCase.mp.path.falco filesep 'models']));
            addpath(genpath([testCase.mp.path.falco filesep 'setup']));
            addpath(genpath([testCase.mp.path.falco filesep 'lib']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco filesep 'models']))
            rmpath(genpath([testCase.mp.path.falco filesep 'setup']))
            rmpath(genpath([testCase.mp.path.falco filesep 'lib']));
        end
    end
    
%% *Tests*
%
% *Creates tests:*
%
% # *testFunctionalMSWC* Input parameters are predefined in the
% ConfigurationFLC.m function which is called by the properties of the
% class and passed in to the test methods. The code in the test methods
% performs the Wavefront Sensing and Control first, then it test verifies
% Iraw, Iest, Iinco, out.complexProjection, out.log10regHist, dm1vp, and thput
% against defined constraints respectively.
%
    methods (Test)     
        function testFunctionalMSWC(testCase)
            mp=testCase.mp;
            %% Step 3: Perform the Wavefront Sensing and Control
            mp.runLabel = 'test_MSWC';            
            [mp, out] = falco_flesh_out_workspace(mp);
            mp = falco_compute_psf_norm_factor(mp);

            % Apply a grid of spots to the input pupil to allow SNWC
            block0 = ones(5, 5);
            block0(3, 3) = 0.7;
            dotGrid = repmat(block0, [51, 51]);
            dotGrid = pad_crop(dotGrid, mp.P1.compact.Narr, 'extrapval', 1);
            % figure; imagesc(dotGrid); axis xy equal tight; colorbar;
            for si = 1:mp.Nsbp
                wvl = mp.sbp_centers(si);
                mp.P1.compact.E(:, :, si) = dotGrid .* ones(mp.P1.compact.Narr);
                for wi = 1:mp.Nwpsbp
                    mp.P1.full.E(:, :, wi, si) = dotGrid .* ones(mp.P1.full.Narr);
                end
            end
            
            [mp, out] = falco_wfsc_loop(mp, out);

            %% Tests:
            Iraw = out.IrawCorrHist(1); % 1.5640e-06
            testCase.verifyGreaterThan(Iraw, 1.55e-6)
            testCase.verifyLessThan(Iraw, 1.57e-6)           
            
            Iest = out.IestScoreHist(1); % 3.4768e-07
            testCase.verifyGreaterThan(Iest, 3.46e-7)
            testCase.verifyLessThan(Iest, 3.49e-7)
            
            Iinco = out.IincoCorrHist(1); % 2.7803e-6
            testCase.verifyGreaterThan(Iinco, 2.77e-6)
            testCase.verifyLessThan(Iinco, 2.79e-6)
            
            dm1pv = out.dm1.Spv(1); % 5.5979e-09
            testCase.verifyGreaterThan(dm1pv, 5.5e-9)
            testCase.verifyLessThan(dm1pv, 5.7e-9)
            
            thput = out.thput(end); % 0.2843
            testCase.verifyGreaterThan(thput, 0.2842)
            testCase.verifyLessThan(thput, 0.2844)
            
            testCase.verifyEqual(out.log10regHist(1), -2)
            
        end       
    end    
end





