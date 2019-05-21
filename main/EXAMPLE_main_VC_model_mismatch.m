% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform a VC simple design run.
%  1) Load the default model parameters for a vortex.
%  2) Specify the values to overwrite.
%  3) Run a single trial of WFC using FALCO.
%
% REVISION HISTORY:
% --------------
% Modified on 2019-04-22 by A.J. Riggs to re-organize the system ID code.
% Modified on 2019-04-12 by He Sun to add the system ID function.
% Modified on 2019-02-26 by A.J. Riggs to load the defaults first.
% ---------------

clear all;
% clc
%% Tell Matlab which Python to use from where. Needed for the E-M algorithm.
% setenv('PATH', '/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin') %--Put the parent directory of the python version you want Matlab to use.
% pyversion('/usr/local/opt/python3/bin/python3.6'); %--Put the location of the python version you want Matlab to use.

%% Step 1: Define Necessary Paths on Your Computer System

% %--Library locations. FALCO and PROPER are required. CVX is optional.
% mp.path.falco = '/Users/ajriggs/Repos/falco-matlab/';  %--Location of FALCO
% mp.path.proper = '/Users/ajriggs/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% % mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX
% 
% %%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = '/Users/ajriggs/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
% mp.path.ws = '/Users/ajriggs/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = '/home/hcst/falco-matlab';%'~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '/home/hcst/PROPER';%'~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '/home/hcst/falco-matlab/data/configs';%'~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '/home/hcst/falco-matlab/data/ws';%'~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path
% rmpath([mp.path.cvx 'lib/narginchk_:']) %--Legend plotting issue if CVX's narginchk function is used instead of Matlab's default function.


%% Step 2: Load default model parameters

EXAMPLE_defaults_VC_simple

%% Step 3: Overwrite default values as desired

%--Record Keeping
mp.SeriesNum = 40;
mp.TrialNum = 2;

%--WFSC Iterations and Control Matrix Relinearization
mp.controller = 'gridsearchEFC';
mp.Nitr = 10; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1;%1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.dm_ind = [1 2]; %--Which DMs to use
mp.ctrl.log10regVec = -1;%--log10 of the regularization exponents (often called Beta values)

%--Training the model
mp.flagTrainModel = true;%false;%
mp.est.flagUseJac = true;%false;%
mp.NitrTrain = mp.Nitr; %--How many iterations to use per training set.

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;
mp.flagUseLearnedJac = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;%k_runTrial;%

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

%--DEBUGGING IN MONOCHROMATIC LIGHT
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation


%% Tuning parameters for System Identification
flagUseTensorflow = true;
if flagUseTensorflow
    mp.est.lr  = 5e-8; % learning rate
    mp.est.lr2 = 1e-2; % learning rate2
    mp.est.epoch = 10; % 
    mp.est.Q0 = 1e-9;
    mp.est.Q1 = 0.4;
    mp.est.R0 = 1e-25;
    mp.est.R1 = 1e-25;
    mp.est.R2 = 0.4;
%     setenv('PATH', '/home/hcst/anaconda3/envs/py35/bin')
%     pyversion('/home/hcst/anaconda3/envs/py35/bin/python3.5')
%     systemID_mod = py.importlib.import_module('falco_systemID');
else
    mp.est.EMitr = 2;
    mp.est.Q = 3e-9;
    mp.est.R = 3e-14;
    mp.est.delta1 = 1e-1;
    mp.est.delta2 = 1e-2;
end
mp.est.flagUseTensorflow = flagUseTensorflow;
%% Sources of model mismatch to include in full model

% %--Generate and save the errors the first time.
% DM1V = 2*randn(34);
% DM2V = 2*randn(34);
% save('/home/hcst/falco-matlab/data/maps/dm_errors.mat','DM1V','DM2V')

%--Load the errors the following times
% load('/Users/ajriggs/Repos/falco-matlab/data/maps/dm_errors.mat','DM1V','DM2V')
load('/home/hcst/falco-matlab/data/maps/dm_errors.mat','DM1V','DM2V')
mp.full.dm1.V0 = DM1V;
mp.full.dm2.V0 = DM2V;

%--Mis-align the DMs for added errors
mp.full.dm1.xc = 16.5 - 0.5; %--16.5 is centered and expected
mp.full.dm2.yc = 16.5 + 0.5; %--16.5 is centered and expected


%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control
% [out, data_train] = falco_adaptive_wfsc_loop(mp);

for klearning = 1 : 10
    if klearning == 1
        mp.flagUseLearnedJac = false;
    else
        mp.flagUseLearnedJac = true;
    end

    out = falco_wfsc_loop(mp);

    if(mp.flagPlot)
        figure(300); semilogy(0:mp.Nitr,out.InormHist,'Linewidth',3); grid on; set(gca,'Fontsize',20);
        hold on
        drawnow
    end
    if klearning == 1
        lr = mp.est.lr;
        lr2 = mp.est.lr2;
        epoch = mp.est.epoch;
        Q0 = mp.est.Q0;
        Q1 = mp.est.Q1;
        R0 = mp.est.R0;
        R1 = mp.est.R1;
        R2 = mp.est.R2;
    else
        lr = mp.est.lr;
        lr2 = mp.est.lr2;
        epoch = mp.est.epoch;
        jacLearned = load(['/home/hcst/falco-matlab/data/jac/', 'jacStructLearned.mat']);
        Q0 = jacLearned.noise_coef(1);
        Q1 = jacLearned.noise_coef(2);
        R0 = jacLearned.noise_coef(3);
        R1 = jacLearned.noise_coef(4);
        R2 = jacLearned.noise_coef(5);
    end
    pycommand = ['python falco_systemID_main.py /home/hcst/falco-matlab/data/jac/ ', ...
                num2str(lr), ' ', num2str(lr2), ' ', num2str(epoch), ' ', 'True ', ...
                num2str(Q0), ' ', num2str(Q1), ' ', num2str(R0), ' ', num2str(R1), ' ', ...
                num2str(R2)];
    
    system(pycommand);
end



