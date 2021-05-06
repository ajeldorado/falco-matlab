%% Functional Test Vortex
%
% The test script will perform the Wavefront Sensing and Control first, then it will
% test verify the values of the outputs to expected values define for each
% test.
classdef TestFunctionalVortex < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we use the mp.path.falco to addpath to the
% function being tested.
    properties
        mp=ConfigurationVortex();
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

    methods (Test)     
        function testFunctionalVortex(testCase)
            mp = testCase.mp;
            
            %% Step 3: Perform the Wavefront Sensing and Control
            mp.runLabel = 'test_vortex';
            [mp, out] = falco_flesh_out_workspace(mp);
            [mp, out] = falco_wfsc_loop(mp, out);

            %% Tests:
            Iend = out.IrawCorrHist(end); % 1.38e-10
            testCase.verifyGreaterThan(Iend, 1.2e-10)
            testCase.verifyLessThan(Iend, 1.5e-10) 
            
            dm1pv = out.dm1.Spv(end); % 1.5057e-08
            testCase.verifyGreaterThan(dm1pv, 1.4e-8)
            testCase.verifyLessThan(dm1pv, 1.6e-8)
            
            thput = out.thput(end); % 28.52%
            testCase.verifyGreaterThan(thput, 0.28)
            testCase.verifyLessThan(thput, 0.29)
            
            import matlab.unittest.constraints.EveryElementOf
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(out.log10regHist, IsEqualTo([-4.5; -4; -4])) 
        end       
    end    
end





