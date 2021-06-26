%% Functional Test FLC
%
% The test script will perform the Wavefront Sensing and Control first, then it will
% test verify the values of the outputs to expected values define for each
% test.
classdef TestFunctionalFLC < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we use the mp.path.falco to addpath to the
% function being tested.
    properties
        mp=ConfigurationFLC();
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath([testCase.mp.path.falco filesep 'setup']));
            addpath(genpath([testCase.mp.path.falco filesep 'lib']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
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
%
    methods (Test)     
        function testFunctionalFLC(testCase)
            mp=testCase.mp;
           
            % Perform the Wavefront Sensing and Control
            mp.runLabel = 'test_FLC';
            [mp, out] = falco_flesh_out_workspace(mp);
            [mp, out] = falco_wfsc_loop(mp, out);
            
            % Tests:
            Iraw = out.IrawCorrHist(end); % 9.9421e-06
            %9.9e-6 < Iraw &&  Iraw < 1e-5
            testCase.verifyGreaterThan(Iraw,9.9e-6)
            testCase.verifyLessThan(Iraw,1e-5)
            
            Iest = out.IestScoreHist(end); % 8.3500e-06
            %8.3e-6 < Iest && Iest < 8.4e-6
            testCase.verifyGreaterThan(Iest,8.3e-6)
            testCase.verifyLessThan(Iest,8.4e-6)
            
            Iinco = out.IincoCorrHist(end); % 1.1000e-05
            %1.0e-5 < Iinco && Iinco < 1.2e-5
            testCase.verifyGreaterThan(Iinco,1.0e-5)
            testCase.verifyLessThan(Iinco,1.2e-2)
            
            complexProj = out.complexProjection(2,1); % 0.7668
            testCase.verifyGreaterThan(complexProj, 0.76)
            testCase.verifyLessThan(complexProj, 0.77)
            
            dm1pv = out.dm1.Spv(end); % 5.8077e-08
            %5.80e-8 < dm1pv && dm1pv < 5.82e-8
            testCase.verifyGreaterThan(dm1pv,5.80e-8)
            testCase.verifyLessThan(dm1pv,5.82e-8)
            
            thput = out.thput(end); % 0.1496
            %0.149 < thput && thput < 0.15
            testCase.verifyGreaterThan(thput,0.149)
            testCase.verifyLessThan(thput,0.15);
            
            %all(out.log10regHist == [-2; -2; -2])
            import matlab.unittest.constraints.EveryElementOf
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(EveryElementOf(out.log10regHist), IsEqualTo(-2))            
        end       
    end    
end





