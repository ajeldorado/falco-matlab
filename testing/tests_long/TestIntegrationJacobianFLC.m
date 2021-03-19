%% Functional Test Jacobian FLC
%
% The test script will perform the Wavefront Sensing and Control first, then 
% it will compute the Jacobians via differencing, and then it will
% test verify the values of the outputs to expected values define for each
% test.
classdef TestIntegrationJacobianFLC < matlab.unittest.TestCase
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
% # *testFunctionalJacobianFLC* Input parameters are predefined in the
% ConfigurationFLC.m function which is called by the properties of the
% class and passed in to the test methods. The code in the test methods
% performs the Wavefront Sensing and Control first, then it test verifies
% rmsNormErrorDM1, and rmsNormErrorDM2 against defined constraints respectively.
%
    methods (Test)     
        function testIntegrationFLC(testCase)
            mp=testCase.mp;
           
            %% Step 3: Perform the Wavefront Sensing and Control
            mp.runLabel = 'test_LC';          
            [mp, out] = falco_flesh_out_workspace(mp);
            mp.dm1.V = zeros(mp.dm1.Nact);
            mp.dm2.V = zeros(mp.dm2.Nact);
            jacStruct = model_Jacobian(mp); %--Get structure containing Jacobians
            
            %%
            G1fastAll = jacStruct.G1;
            G2fastAll = jacStruct.G2;
            Nind = 20;
            subinds = 4*(1:Nind);
            absG1sum = sum(abs(G1fastAll), 1);
            indG1 = find(absG1sum > 1e-2*max(absG1sum));
            indG1subset = indG1(subinds); % Take a 20-actuator subset
            absG2sum = sum(abs(G2fastAll), 1);
            indG2 = find(absG2sum > 1e-2*max(absG2sum));
            indG2subset = indG2(subinds); % Take a 20-actuator subset
            G1fast = G1fastAll(:, indG1subset);
            G2fast = G2fastAll(:, indG2subset);
            
            %% Compute Jacobian via differencing (slower)
            modvar.whichSource = 'star';
            modvar.sbpIndex = 1;
            modvar.starIndex = 1;
            Eunpoked2D = model_compact(mp, modvar);
            Eunpoked = Eunpoked2D(mp.Fend.corr.maskBool);
            DeltaV = 1e-4;
            G1slow = zeros(mp.Fend.corr.Npix, Nind);
            G2slow = zeros(mp.Fend.corr.Npix, Nind);
            
            for ii = 1:Nind
                % DM1
                mp.dm1.V = zeros(mp.dm1.Nact);
                mp.dm2.V = zeros(mp.dm2.Nact);
                mp.dm1.V(indG1subset(ii)) = DeltaV;
                Epoked2D = model_compact(mp, modvar);
                Epoked = Epoked2D(mp.Fend.corr.maskBool);
                G1slow(:, ii) = (Epoked - Eunpoked) / DeltaV;
                
                % DM2
                mp.dm1.V = zeros(mp.dm1.Nact);
                mp.dm2.V = zeros(mp.dm2.Nact);
                mp.dm2.V(indG2subset(ii)) = DeltaV;
                Epoked2D = model_compact(mp, modvar);
                Epoked = Epoked2D(mp.Fend.corr.maskBool);
                G2slow(:, ii) = (Epoked - Eunpoked) / DeltaV; 
            end
            
            %% Tests
            
            rmsNormErrorDM1=sqrt(sum(abs(G1slow(:) - G1fast(:)).^2)/sum(abs(G1slow(:)).^2));
            rmsNormErrorDM2 = sqrt(sum(abs(G2slow(:) - G2fast(:)).^2)/sum(abs(G2slow(:)).^2));
            testCase.verifyLessThan(rmsNormErrorDM1,1e-3)
            testCase.verifyLessThan(rmsNormErrorDM2,1e-3)          
        end       
    end    
end





