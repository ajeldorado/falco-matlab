%% Test DM Constraint Functions
%
% Unit tests of the methods in the ConstrainDM class.
%
classdef TestActLimits < matlab.unittest.TestCase
%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    properties
%         increasing_data = reshape(0:24, [5, 5]).';
%         data_to_fail = ones(5) + eye(5);
%         data_ones = ones(5);
%         data_end = zeros(12, 12);
%         bigeye = 5.0 * eye(5);
%         tie_no_tied_dead = zeros(3, 3);
%         volts_no_tied_dead = reshape(1:9, [3, 3]).';
%         tie_tied = [0, 1, 0; 0, 1, 3; 2, 2, 0];
%         volts_tied = [1, 2, 3; 4, 2, 6; 7, 7, 9];
%         tie_dead = [-1, 0, -1; 0, 0, 0; 0, 0, 0];
%         volts_dead = [0, 2, 0; 4, 5, 6; 7, 8, 9];
%         tie_tied_dead = [-1, 1, -1; 0, 1, 3; 2, 2, 0];
%         volts_tied_dead = [0, 2, 0; 4, 2, 6; 7, 7, 9];
%         
%         rng_seed = 2021;
%         nact = 48;
%         flatmap_flat = zeros(48);
%         vmin = 0.;
%         vmax = 100.;
%         vmid = 50
%         vlat = 50.;
%         vdiag = 75.;
%         vquant = 100 / 2^16;
%         maxiter = 1000;
%         
%         vneighbor = 30.;
%         vcorner = 30.;
%         dm0 = 50*ones(48);
%         uniform_flat = 50*ones(48);
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
    
%% Unit tests of ConstrainDM

    methods (Test)
                
        function test_grow(testCase)
            % Test tied-group growth by adding in a chain of NR actuators which all
            % connect along different axes but should be linked together by the
            % grow_group() loop
            
            dmvobj = DM;
            
            Vmin = 0; %self.dmvobj.Vmin
            Vmax = 100; %self.dmvobj.Vmax
            Nact = 48; %cfg.dmlist[0].registration['nact']

            dmvobj.Vmin = 0;
            dmvobj.Vmax = 100;
            dmvobj.dVnbrLat = (Vmax - Vmin)/4.;
            dmvobj.dVnbrDiag = (Vmax - Vmin)/4.;
            dmvobj.facesheetFlatmap = zeros(Nact);
            dmvobj.VtoH = 4e-9*ones(Nact);
            dmvobj.marginNbrRule = 0.001678466796875;

            % keep these small so we don't accidentally freeze
            % self.dmvobj.dVnbrLat = (Vmax - Vmin)/4.
            % dVnbrLat = self.dmvobj.dVnbrLat
            dVnbrLat = (Vmax - Vmin)/4.;
            % self.dmvobj.dVnbrDiag = (Vmax - Vmin)/4.
            % dVnbrDiag = self.dmvobj.dVnbrDiag
            dVnbrDiag = (Vmax - Vmin)/4.;
            vmid = (Vmin + Vmax)/2.; % midpoint

            dmtest = vmid * ones(Nact, Nact);
            dmtest(1, 1) = vmid + dVnbrLat/2.;
            dmtest(1, 2) = vmid - dVnbrLat/2.;
            dmtest(2, 3) = vmid - dVnbrLat/2. + dVnbrDiag;
            dmtest(3, 4) = vmid - dVnbrLat/2.;
            dmtest(4, 4) = vmid + dVnbrLat/2.;
            dmtest(5, 3) = vmid + dVnbrLat/2. - dVnbrDiag;
            dmtest(6, 2) = vmid + dVnbrLat/2.;
            dmtest(6, 1) = vmid - dVnbrLat/2.;

            tiemapExpected = zeros(size(dmtest));
            tiemapExpected(dmtest ~= 50) = 1;
            
%             figure(1); imagesc(dmtest); axis xy equal tight; colorbar;
            % figure(2); imagesc(B); axis xy equal tight; colorbar;
            % B = reshape(1:nact^2, [nact, nact]);

            [freezeVec, linkList] = ActLimits.maplimits(dmtest, dmvobj);

            tiemap = zeros(size(dmtest));
            for ii = 1:length(linkList)
                tiemap(linkList{ii}) = ii;
            end
%             figure(3); imagesc(tiemap); axis xy equal tight; colorbar;
%             figure(4); imagesc(tiemapExpected); axis xy equal tight; colorbar;

            testCase.verifyTrue(all(all(tiemap == tiemapExpected)))
            testCase.verifyTrue(isempty(freezeVec))

        end

        
        function test_tie_multiple(testCase)
            % Test tied-group growth by adding in a chain of NR actuators which all
            % connect along different axes but should be linked together by the
            % grow_group() loop
            
            dmvobj = DM;
            
            Vmin = 0; %self.dmvobj.Vmin
            Vmax = 100; %self.dmvobj.Vmax
            Nact = 48; %cfg.dmlist[0].registration['nact']

            dmvobj.Vmin = 0;
            dmvobj.Vmax = 100;
            dmvobj.dVnbrLat = 50; %(Vmax - Vmin)/4.;
            dmvobj.dVnbrDiag = 75; %(Vmax - Vmin)/4.;
            dmvobj.facesheetFlatmap = zeros(Nact);
            dmvobj.VtoH = 4e-9*ones(Nact);
            dmvobj.marginNbrRule = 0.001678466796875;

            vneighbor = dmvobj.dVnbrDiag; %(Vmax - Vmin)/4.;
            vcorner = dmvobj.dVnbrDiag; %(Vmax - Vmin)/4.;
            vmid = (Vmin + Vmax)/2.; % midpoint

            dmtest = vmid * ones(Nact, Nact);
            dmtest(2, 2) = vmid + vcorner/2.;
            dmtest(3, 1) = vmid - vcorner/2.;
            dmtest(6, 2) = vmid + vcorner/2.;
            dmtest(7, 3) = vmid - vcorner/2.;
            dmtest(2, 6) = vmid + vneighbor/2.;
            dmtest(3, 6) = vmid - vneighbor/2.;
            dmtest(5, 6) = vmid + vneighbor/2.;
            dmtest(6, 7) = vmid - vneighbor/2.;

            tiemapExpected = zeros(size(dmtest));
            tiemapExpected(2, 2) = 4;
            tiemapExpected(3, 1) = 4;
            tiemapExpected(6, 2) = 2;
            tiemapExpected(7, 3) = 2;
            tiemapExpected(2, 6) = 1;
            tiemapExpected(3, 6) = 1;
            tiemapExpected(5, 6) = 3;
            tiemapExpected(6, 7) = 3;
            
%             figure(1); imagesc(dmtest); axis xy equal tight; colorbar;

            [freezeVec, linkList] = ActLimits.maplimits(dmtest, dmvobj);

            tiemap = zeros(size(dmtest));
            for ii = 1:length(linkList)
                tiemap(linkList{ii}) = ii;
            end
%             figure(3); imagesc(tiemap); axis xy equal tight; colorbar;
%             figure(4); imagesc(tiemapExpected); axis xy equal tight; colorbar;

            testCase.verifyTrue(all(all(tiemap == tiemapExpected)))
            testCase.verifyTrue(isempty(freezeVec))
        
        end

        
        function test_tie_dead_grow(testCase)
            % Test tied-group growth by adding in a chain of NR actuators which all
            % connect along different axes but should be linked together by the
            % grow_group() loop
            
            dmvobj = DM;
            
            Vmin = 0; %self.dmvobj.Vmin
            Vmax = 100; %self.dmvobj.Vmax
            Nact = 48; %cfg.dmlist[0].registration['nact']

            dmvobj.Vmin = 0;
            dmvobj.Vmax = 100;
            dmvobj.dVnbrLat = 50; %(Vmax - Vmin)/4.;
            dmvobj.dVnbrDiag = 75; %(Vmax - Vmin)/4.;
            dmvobj.facesheetFlatmap = zeros(Nact);
            dmvobj.VtoH = 4e-9*ones(Nact);
            dmvobj.marginNbrRule = 0.001678466796875;

            vneighbor = dmvobj.dVnbrDiag; %(Vmax - Vmin)/4.;
            vcorner = dmvobj.dVnbrDiag; %(Vmax - Vmin)/4.;
            vmid = (Vmin + Vmax)/2.; % midpoint
            
            tmptie = zeros(Nact, Nact);
            tmptie(2, 7) = -1;

            dmtest = vmid * ones(Nact, Nact);
            dmtest(1, 1) = vmid + vneighbor/2.;
            dmtest(1, 2) = vmid - vneighbor/2.;
            dmtest(2, 2) = vmid + vneighbor/2.; % 0,0 & 1,1 same, so no diag viol
            dmtest(2, 3) = vmid - vneighbor/2.;
            dmtest(2, 4) = vmid + vneighbor/2.;
            dmtest(2, 5) = vmid - vneighbor/2.;
            dmtest(2, 6) = vmid + vneighbor/2.;
            dmtest(2, 7) = vmid - vneighbor/2.;
            
            freezemapExpected = zeros(size(dmtest));
            freezemapExpected(1, 1) = -1;
            freezemapExpected(1, 2) = -1;
            freezemapExpected(2, 2) = -1;
            freezemapExpected(2, 3) = -1;
            freezemapExpected(2, 4) = -1;
            freezemapExpected(2, 5) = -1;
            freezemapExpected(2, 6) = -1;
            freezemapExpected(2, 7) = -1;
            
%             figure(1); imagesc(dmtest); axis xy equal tight; colorbar;

            [freezeVec, linkList] = ActLimits.maplimits(dmtest, dmvobj, tmptie);
            disp(linkList)
            disp(length(linkList))

            freezemap = zeros(size(dmtest));
            freezemap(freezeVec) = -1;

%             figure(3); imagesc(freezemap); axis xy equal tight; colorbar;
%             figure(4); imagesc(freezemapExpected); axis xy equal tight; colorbar;

            testCase.verifyTrue(all(all(freezemap == freezemapExpected)))
            testCase.verifyTrue(isempty(linkList))
        
        end        
        
    end
    
end
