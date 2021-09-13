%% Integration Test LC
%
classdef TestIntegrationJacobianLC < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters is loaded to be used by the methods.
%
    properties
        mp = ConfigurationLC();
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
        
        function testJacobianLC(testCase)
            mp = testCase.mp;
            mp.runLabel = 'test_LC';
            [mp, out] = falco_flesh_out_workspace(mp);
            
            %% Fast Jacobian calculation
            mp.dm1.V = zeros(mp.dm1.Nact);
            mp.dm2.V = zeros(mp.dm2.Nact);
            jacStruct = model_Jacobian(mp); %--Get structure containing Jacobians
            
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
            rmsNormErrorDM1 = sqrt(sum(abs(G1slow(:) - G1fast(:)).^2)/sum(abs(G1slow(:)).^2));
            rmsNormErrorDM2 = sqrt(sum(abs(G2slow(:) - G2fast(:)).^2)/sum(abs(G2slow(:)).^2));
            
            testCase.verifyLessThan(rmsNormErrorDM1,1e-3)
            testCase.verifyLessThan(rmsNormErrorDM2,1e-3)
            
        end
        
        function testJacobianLCmftAS(testCase)
            mp = testCase.mp;
            mp.runLabel = 'test_LC';
            mp.propMethodPTP = 'mft';
            [mp, out] = falco_flesh_out_workspace(mp);
            
            %% Fast Jacobian calculation
            mp.dm1.V = zeros(mp.dm1.Nact);
            mp.dm2.V = zeros(mp.dm2.Nact);
            jacStruct = model_Jacobian(mp); %--Get structure containing Jacobians
            
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
            rmsNormErrorDM1 = sqrt(sum(abs(G1slow(:) - G1fast(:)).^2)/sum(abs(G1slow(:)).^2));
            rmsNormErrorDM2 = sqrt(sum(abs(G2slow(:) - G2fast(:)).^2)/sum(abs(G2slow(:)).^2));
            
            testCase.verifyLessThan(rmsNormErrorDM1, 0.02)
            testCase.verifyLessThan(rmsNormErrorDM2, 0.02)
            
        end
        
        function testJacobianLCnoFPM(testCase)
            mp = testCase.mp;
            mp.jac.minimizeNI = true;
            
            mp.runLabel = 'test_LC';
            [mp, out] = falco_flesh_out_workspace(mp);
            mp.dm1.V = zeros(mp.dm1.Nact);
            mp.dm2.V = zeros(mp.dm2.Nact);

            %% Compute Jacobian with the Jacobian model
            %--Calculate the starting DM surfaces beforehand to save time.
            if(any(mp.dm_ind==1)); mp.dm1.compact.surfM = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.NdmPad); else; mp.dm1.compact.surfM = zeros(2); end
            if(any(mp.dm_ind==2)); mp.dm2.compact.surfM = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.NdmPad); else; mp.dm2.compact.surfM = zeros(2); end

            iMode = 1;
            G1fastAll = model_Jacobian_no_FPM(mp, iMode, 1);
            G2fastAll = model_Jacobian_no_FPM(mp, iMode, 2);

            Nind = 20;
            subinds = 4*(1:Nind);
            absG1sum = sum(abs(G1fastAll), 1);
            thresh = 1e-1;
            indG1 = find(absG1sum > thresh*max(absG1sum));
            indG1subset = indG1(subinds); % Take a 20-actuator subset
            absG2sum = sum(abs(G2fastAll), 1);
            indG2 = find(absG2sum > thresh*max(absG2sum));
            indG2subset = indG2(subinds); % Take a 20-actuator subset
            G1fast = G1fastAll(:, indG1subset);
            G2fast = G2fastAll(:, indG2subset);

            %% Compute Jacobian via differencing with the compact model (slower)

            % Get the unocculted peak E-field and coronagraphic E-field
            if mp.jac.minimizeNI
                modvar.sbpIndex = mp.jac.sbp_inds(iMode);
                modvar.zernIndex = mp.jac.zern_inds(iMode);
                modvar.starIndex = mp.jac.star_inds(iMode);
                modvar.whichSource = 'star';
                Eunocculted = model_compact(mp, modvar, 'nofpm');
                [~, indPeak] = max(abs(Eunocculted(:)));
            end

            modvar.whichSource = 'star';
            modvar.sbpIndex = 1;
            modvar.starIndex = 1;
            Eunpoked2D = model_compact(mp, modvar, 'nofpm');
            Eunpoked = Eunpoked2D(indPeak);
            DeltaV = 1e-4;
            G1slow = zeros(1, Nind);
            G2slow = zeros(1, Nind);

            for ii = 1:Nind
                % DM1
                mp.dm1.V = zeros(mp.dm1.Nact);
                mp.dm2.V = zeros(mp.dm2.Nact);
                mp.dm1.V(indG1subset(ii)) = DeltaV;
                Epoked2D = model_compact(mp, modvar, 'nofpm');
                Epoked = Epoked2D(indPeak);
                G1slow(:, ii) = (Epoked - Eunpoked) / DeltaV;

                % DM2
                mp.dm1.V = zeros(mp.dm1.Nact);
                mp.dm2.V = zeros(mp.dm2.Nact);
                mp.dm2.V(indG2subset(ii)) = DeltaV;
                Epoked2D = model_compact(mp, modvar, 'nofpm');
                Epoked = Epoked2D(indPeak);
                G2slow(:, ii) = (Epoked - Eunpoked) / DeltaV;   
            end

            %% Tests
            rmsNormErrorDM1 = sqrt(sum(abs(G1slow(:) - G1fast(:)).^2)/sum(abs(G1slow(:)).^2));
            rmsNormErrorDM2 = sqrt(sum(abs(G2slow(:) - G2fast(:)).^2)/sum(abs(G2slow(:)).^2));
            
            testCase.verifyLessThan(rmsNormErrorDM1,1e-3)
            testCase.verifyLessThan(rmsNormErrorDM2,1e-3)
        end
        
    end    
end





