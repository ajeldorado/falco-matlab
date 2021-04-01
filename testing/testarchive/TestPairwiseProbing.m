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
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
    properties
        mp=Parameters();
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath([testCase.mp.path.falco filesep 'lib']));
            addpath(genpath([testCase.mp.path.falco filesep 'models']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco filesep 'lib']));
            rmpath(genpath([testCase.mp.path.falco filesep 'models']));
        end
    end
    
%% Tests
%
%  Creates four tests:
%
% # *testPairwiseProbing* Test falco_est_pairwise_probing.m.                                   
    methods (Test)    
        function testPairwiseProbing(testCase)
            
            EXAMPLE_defaults_try_running_FALCO
            
            mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
            mp.flagPlot = false;
            mp.SeriesNum = 1;
            mp.TrialNum = 1;
            mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
            mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
            mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
            mp.Nitr = 3; %--Number of wavefront control iterations
            mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
                mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
                '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
                '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
                '_',mp.controller];
            
            mp.whichPupil = 'SIMPLE';
            mp.P1.ODnorm = 1;
            [mp, out] = falco_flesh_out_workspace(mp);
            N = size(mp.P1.full.E, 1);
            alpha = 2.5;
            mirror_figure = 1e-9;
            errormap = falco_gen_simple_PSD_errormap(N, alpha, mirror_figure);
            mp.P1.full.E = exp(2*pi*1j/mp.lambda0*errormap);
            mp.P1.compact.E = exp(2*pi*1j/mp.lambda0*errormap);
            Im = falco_get_summed_image(mp);
            mp = falco_compute_PSF_norm_factor(mp);
            mp.estimator = 'perfect';
            ev = falco_est_perfect_Efield_with_Zernikes(mp);
            Etrue = ev.Eest;
            mp.estimator = 'pwp-bp-square';
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
        end
    end    
end