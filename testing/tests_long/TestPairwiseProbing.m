%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test Pairwise Probing
%
% We test the function falco_est_pairwise_probing.m.
classdef TestPairwiseProbing < matlab.unittest.TestCase    
%% Properties
%
%% Properties
%
% A presaved file with FALCO parameters was saved and is loaded to be used
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
            addpath(genpath('../../config'));
            addpath(genpath('../../lib'));
            addpath(genpath('../../lib_external'));
            addpath(genpath('../../models'));
            addpath(genpath('../../setup'));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath('../../config'));
            rmpath(genpath('../../lib'));
            rmpath(genpath('../../lib_external'));
            rmpath(genpath('../../models'));
            rmpath(genpath('../../setup'));
        end
    end
    
%% Tests
%
%  Creates four tests:
%
% # *testPairwiseProbing* Test falco_est_pairwise_probing.m.                                   
    methods (Test)    
        function testPairwiseProbing(testCase)
            
%             EXAMPLE_defaults_try_running_FALCO
            mp = testCase.mp;
            
            mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
            mp.flagPlot = false;
            mp.SeriesNum = 1;
            mp.TrialNum = 1;
            mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
            mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
            mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
            mp.Nitr = 3; %--Number of wavefront control iterations
            mp.runLabel = 'testing_probing';
            
            % Overwrite the LUVOIR-B pupil with an open circle to get
            % better contrast for this test.
            % Inputs common to both the compact and full models
            inputs.ID = 0;
            inputs.OD = 1.0;
            % Full model
            inputs.Nbeam = mp.P1.full.Nbeam;
            inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
            mp.P1.full.mask = falco_gen_pupil_Simple(inputs); 
            % Compact model
            inputs.Nbeam = mp.P1.compact.Nbeam;
            inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
            mp.P1.compact.mask = falco_gen_pupil_Simple(inputs); 

            [mp, out] = falco_flesh_out_workspace(mp);
            N = size(mp.P1.full.E, 1);
            alpha = 2.5;
            mirror_figure = 1e-10;
            errormap = falco_gen_simple_PSD_errormap(N, alpha, mirror_figure);
            mp.P1.full.E = exp(2*pi*1j/mp.lambda0*errormap);
            mp.P1.compact.E = exp(2*pi*1j/mp.lambda0*errormap);
            Im = falco_get_summed_image(mp);
            mp = falco_compute_psf_norm_factor(mp);
            
            % Get exact E-field for comparison:
            mp.estimator = 'perfect';
            ev = falco_est_perfect_Efield_with_Zernikes(mp);
            Etrue = ev.Eest;
            
            % Estimate E-field with square-defined probing region
            mp.estimator = 'pairwise';
            mp.est = rmfield(mp.est, 'probe');
            mp.est.probe = Probe;
            mp.est.probe.Npairs = 3;     % Number of pair-wise probe PAIRS to use.
            mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
            mp.est.probe.radius = 12;    % Max x/y extent of probed region [lambda/D].
            mp.est.probe.xOffset = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
            mp.est.probe.yOffset = 0;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
            mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
            mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.
            ev.dummy = 1;
            ev = falco_est_pairwise_probing(mp, ev);
            Eest = ev.Eest;
            Eest2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
            Eest2D(mp.Fend.corr.maskBool) = Eest;
            Etrue2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
            Etrue2D(mp.Fend.corr.maskBool) = Etrue;
            Etrue2D(Eest2D == 0) = 0;
            meanI = mean(abs(Etrue).^2);
            meanIdiff = mean(abs(Etrue2D(mp.Fend.corr.maskBool)-Eest2D(mp.Fend.corr.maskBool)).^2);
            percentEstError = meanIdiff/meanI*100;
            testCase.verifyLessThan(percentEstError, 4.0)
            
            % Estimate E-field with rectangle-defined probing region
            mp.estimator = 'pairwise-rect';
            mp.est = rmfield(mp.est, 'probe');
            mp.est.probe = Probe;
            mp.est.probe.Npairs = 3;     % Number of pair-wise probe PAIRS to use.
            mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
            % mp.est.probe.radius = 12;    % Max x/y extent of probed region [lambda/D].
            mp.est.probe.xOffset = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
            mp.est.probe.yOffset = 0;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
            % mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
            mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.
            mp.est.probe.xiOffset = 6;
            mp.est.probe.etaOffset = 0;
            mp.est.probe.width = 12;        
            mp.est.probe.height = 24;   
            ev.dummy = 1;
            ev = falco_est_pairwise_probing(mp, ev);
            Eest = ev.Eest;
            Eest2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
            Eest2D(mp.Fend.corr.maskBool) = Eest;
            Etrue2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
            Etrue2D(mp.Fend.corr.maskBool) = Etrue;
            Etrue2D(Eest2D == 0) = 0;
            % Block out the unprobed strip in the middle.
            xMiddle = ceil((mp.Fend.Neta+1)/2);
            Eest2D(:, xMiddle-1:xMiddle+1) = 0;
            Etrue2D(:, xMiddle-1:xMiddle+1) = 0;
            meanI = mean(abs(Etrue).^2);
            meanIdiff = mean(abs(Etrue2D(mp.Fend.corr.maskBool)-Eest2D(mp.Fend.corr.maskBool)).^2);
            percentEstError = meanIdiff/meanI*100;
            testCase.verifyLessThan(percentEstError, 4.0)
        end
    end
    
end
