%% Functional Test FLC
%
% The test script will perform the Wavefront Sensing and Control first, then it will
% test verify the values of the outputs to expected values define for each
% test.
classdef TestFunctionalFLC < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is loaded to be used
% by methods. In this case we use the mp.path.falco to addpath to the
% function being tested.
    properties
        mp = ConfigurationFLC();
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
% # *testFunctionalFLC* Input parameters are predefined in the
% ConfigurationFLC.m function which is called by the properties of the
% class and passed in to the test methods. The code in the test methods
% performs the Wavefront Sensing and Control first, then it test verifies
% Iraw, Iest, Iinco, out.complexProjection, dm1vp, thput, and
% out.log10regHist against defined constraints respectively.

    methods (Test)     
        function testFunctionalFLC(testCase)
            mp=testCase.mp;
           
            % Perform the Wavefront Sensing and Control
            mp.runLabel = 'test_FLC_WFSC';
            [mp, out] = falco_flesh_out_workspace(mp);
            [mp, out] = falco_wfsc_loop(mp, out);
            
            % Tests:
            Iraw = out.IrawCorrHist(end);
            testCase.verifyTrue(isalmost(Iraw, 8.4796e-06, 3e-7))
            
            Iest = out.IestScoreHist(end);
            testCase.verifyTrue(isalmost(Iest, 2.7282e-05, 3e-6))
            
            Iinco = out.IincoCorrHist(end);
            testCase.verifyTrue(isalmost(Iinco, -1.2982e-05, 3e-6))
            
            complexProj = out.complexProjection(2, 1);
            testCase.verifyTrue(isalmost(complexProj, 1.2403, 1e-2))
            
            dm1pv = out.dm1.Spv(end);
            testCase.verifyTrue(isalmost(dm1pv, 7.7374e-08, 1e-9))

            thput = out.thput(end);
            testCase.verifyTrue(isalmost(thput, 0.1455, 1e-3))
            
            import matlab.unittest.constraints.EveryElementOf
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(EveryElementOf(out.log10regHist), IsEqualTo(-2))            
        end       
    end    
end
