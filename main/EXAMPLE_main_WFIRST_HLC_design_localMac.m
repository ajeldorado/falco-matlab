% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform an HLC design run.
%  1) Load the default model parameters for an HLC.
%  2) Specify the values to overwrite.
%  3) Run a single trial of WFC using FALCO.
%
% REVISION HISTORY:
% --------------
% Modified on 2019-02-26 by A.J. Riggs to load the defaults first.
% ---------------

clear all;
% close all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = '/Users/jllopsay/Documents/GitHub/falco-matlab/';  %--Location of FALCO
mp.path.proper = '/Users/jllopsay/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '/Users/jllopsay/Documents/GitHub/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '/Users/jllopsay/Documents/GitHub/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path
% rmpath([mp.path.cvx 'lib/narginchk_:']) %--Legend plotting issue if CVX's narginchk function is used instead of Matlab's default function.


%% Step 2: Load default model parameters

% EXAMPLE_defaults_WFIRST_HLC_design_local
EXAMPLE_defaults_WFIRST_HLC_design_slowpoke
mp.aux.omega = 0;
mp.aux.gamma = 1e-4;
mp.aux.firstOmegaItr = 11111110;
mp.aux.betaMinusOne = 0;
mp.aux.ItrDump = 0;
mp.aux.flagOmega=0;
mp.aux.omegaMin = 2;
mp.aux.omegaMax = 9;
mp.aux.minNIprob = 1111111;
mp.aux.firstDM9It = 0;
mp.aux.NumConvCFIt = 2;
mp.aux.wDM9_arr = [];
mp.aux.flagRegDM9 = false;
mp.aux.betadm9Min = -8;
mp.aux.betadm9Max = -5;
mp.aux.firstRegDM9Itr = 111110;
SetA2 = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
SetB2 = [1, 1j, 129, 1, 1];

mp.ctrl.sched_mat = [...
     repmat(SetA2,[min([mp.Nitr,mp.aux.firstDM9It]),1]);...
     repmat(SetB2,[mp.Nitr-mp.aux.firstDM9It,1])...
     ];

[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 2;
mp.TrialNum = 5;

%--Force DM9 to be mirror symmetric about y-axis
NactTotal = ceil_even(mp.dm9.actres*mp.F3.Rin*2)^2; %-NOTE: This will be different if influence function for DM9 is not '3x3'. Needs to be the same value as mp.dm9.NactTotal, which is calculated later.
Nact = sqrt(NactTotal);
LinIndMat = zeros(Nact); %--Matrix of the linear indices
LinIndMat(:) = 1:NactTotal;
FlippedLinIndMat = fliplr(LinIndMat);
% mp.dm9.tied = zeros(NactTotal/2,2);
% for jj=1:NactTotal/2
%     mp.dm9.tied(jj,1) = LinIndMat(jj);
%     mp.dm9.tied(jj,2) = FlippedLinIndMat(jj);
% end

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

% %--DEBUGGING ONLY: Monochromatic light
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

%%
mp.dm9.actres = 5;%dm9_actres_arr(randperm(length(dm9_actres_arr),1));
% mp.min_azimSize_dm9 = 5;%min_azimSize_dm9_arr(randperm(length(min_azimSize_dm9_arr),1));
mp.max_azimSize_dm9 = 1;

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);

