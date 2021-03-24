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
            addpath(genpath([testCase.mp.path.falco filesep 'setup']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco filesep 'setup']))
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

            % Apply a grid of spots to the input pupil to allow SNWC
            block0 = ones(5, 5);
            block0(3, 3) = 0.7;%0;
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
            Iraw = out.IrawCorrHist(end); % 7.7579e-07
            %7.7e-7 < Iraw &&  Iraw < 7.8e-7
            testCase.verifyGreaterThan(Iraw,7.7e-7)
            testCase.verifyLessThan(Iraw,7.8e-7)           
            
            Iest = out.IestScoreHist(end); % 9.1479e-08
            %9.1e-8 < Iest && Iest < 9.2e-8
            testCase.verifyGreaterThan(Iest,9.1e-8)
            testCase.verifyLessThan(Iest,9.2e-8)
            
            Iinco = out.IincoCorrHist(end); % 1.7955e-06
            %1.75e-6 < Iinco && Iinco < 1.8e-6
            testCase.verifyGreaterThan(Iinco,1.75e-6)
            testCase.verifyLessThan(Iinco,1.8e-6)

            %all(out.complexProjection(:) > 0.88)
            import matlab.unittest.constraints.EveryElementOf
            import matlab.unittest.constraints.IsGreaterThan
            testCase.verifyThat(EveryElementOf(out.complexProjection(:)), IsGreaterThan(0.88))
            
            dm1pv = out.dm1.Spv(end); % 1.5482e-08
            %1.5e-8 < dm1pv && dm1pv < 1.6e-8
            testCase.verifyGreaterThan(dm1pv,1.5e-8)
            testCase.verifyLessThan(dm1pv,1.6e-8)
            
            thput = out.thput(end); % 0.2837
            %0.28 < thput && thput < 0.29
            testCase.verifyGreaterThan(thput,0.28)
            testCase.verifyLessThan(thput,0.29)
            
            %all(out.log10regHist == [-2; -2; -2])
            import matlab.unittest.constraints.EveryElementOf
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(EveryElementOf(out.log10regHist), IsEqualTo([-2])) 
        end       
    end    
end





