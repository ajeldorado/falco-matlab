%% Functional Test LC
%
clear

mp=ConfigurationLC();

mp.flagPlot = true;

% Step 3: Perform the Wavefront Sensing and Control
mp.runLabel = 'test_LC';
[mp, out] = falco_flesh_out_workspace(mp);
[mp, out] = falco_wfsc_loop(mp, out);

% The test script will perform the Wavefront Sensing and Control first, then it will
% test verify the values of the outputs to expected values define for each
% test.
% classdef TestFunctionalLC < matlab.unittest.TestCase
% %% Properties
% %
% % A presaved file with FALCO parameters was saved and is lodaded to be used
% % by methods. In this case we use the mp.path.falco to addpath to the
% % function being tested.
%     properties
%         mp=ConfigurationLC();
%     end
% 
% %% Setup and Teardown Methods
% %
% %  Add and remove path to utils functions to be tested.
% %
%     methods (TestClassSetup)
%         function addPath(testCase)
%             addpath(genpath([testCase.mp.path.falco filesep 'setup']));
%             addpath(genpath([testCase.mp.path.falco filesep 'lib']));
%         end
%     end
%     methods (TestClassTeardown)
%         function removePath(testCase)
%             rmpath(genpath([testCase.mp.path.falco filesep 'setup']))
%             rmpath(genpath([testCase.mp.path.falco filesep 'lib']));
%         end
%     end    
%     
% %% *Tests*
% %
% % *Creates tests:*
% %
% % # *testFunctionalLC* Input parameters are predefined in the
% % ConfigurationFLC.m function which is called by the properties of the
% % class and passed in to the test methods. The code in the test methods
% % performs the Wavefront Sensing and Control first, then it test verifies
% % Iraw, Iest, out.IincoCorrHist, out.log10regHist, dm1vp, and thput
% % against defined constraints respectively.
% %
%     methods (Test)     
%         function testFunctionalLC(testCase)
%             mp=testCase.mp;
%             %% Step 3: Perform the Wavefront Sensing and Control
%             mp.runLabel = 'test_LC';
%             [mp, out] = falco_flesh_out_workspace(mp);
%             [mp, out] = falco_wfsc_loop(mp, out);
%             
%             % Tests:
%             Iraw = out.IrawCorrHist(end); % 1.3233e-07
%             %1.3e-07 < Iraw &&  Iraw < 1.4e-07
%             testCase.verifyGreaterThan(Iraw,1.3e-07)
%             testCase.verifyLessThan(Iraw,1.4e-07)           
%             
%             Iest = out.IestScoreHist(end); % 8.8397e-07
%             %8.8e-7 < Iest && Iest < 8.9e-7
%             testCase.verifyGreaterThan(Iest,8.8e-7)
%             testCase.verifyLessThan(Iest,8.9e-7)
%             
%             %all(out.IincoCorrHist == 0)
%             import matlab.unittest.constraints.EveryElementOf
%             import matlab.unittest.constraints.IsEqualTo
%             testCase.verifyThat(EveryElementOf(out.IincoCorrHist), IsEqualTo(0)) 
%             
%             dm1pv = out.dm1.Spv(end); % 8.2674e-08
%             %8.2e-08 < dm1pv && dm1pv < 8.3e-08
%             testCase.verifyGreaterThan(dm1pv,8.2e-08)
%             testCase.verifyLessThan(dm1pv,8.3e-08)
%                      
%             thput = out.thput(end); % 0.0740
%             %0.073 < thput && thput < 0.075
%             testCase.verifyGreaterThan(thput,0.073)
%             testCase.verifyLessThan(thput,0.075);
%             
%             %all(out.log10regHist == [-2; -2; -3])
%             import matlab.unittest.constraints.EveryElementOf
%             import matlab.unittest.constraints.IsEqualTo
%             testCase.verifyThat(out.log10regHist, IsEqualTo([-2; -2; -3]))             
%         end       
%     end    
% end





