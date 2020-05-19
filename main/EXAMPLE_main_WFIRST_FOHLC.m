% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform an first-order approximation HLC (FOHLC) design run.

clear all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command:  addpath(full_path_to_proper); savepath;

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false


%% Step 2: Load default model parameters

EXAMPLE_defaults_WFIRST_FOHLC


%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = true;%false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 35;%1;
mp.TrialNum = 4;%1;

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

%--Use just 1 wavelength for initial debugging of code
mp.fracBW = 0.10;%0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 5;%1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control

% mp.F3.Rin = 5;%2.7;%4;%2.7; % maximum radius of inner part of the focal plane mask, in lambda0/D
% mp.F3.RinA = 2.7;%mp.F3.Rin; % inner hard-edge radius of the focal plane mask (lambda0/D). Needs to be <= mp.F3.Rin 

%--Use just 1 wavelength for initial debugging of code
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation


%--Zernikes to suppress with controller
mp.jac.zerns = 1;%[1,2,3];  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = [1, 1e-9, 1e-9];%1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
  

mp.dm8.weight = 1e-2;%1e1;%1; % Jacobian weight for the FPM dielectric. Smaller weight makes stroke larger by the inverse of this factor.
mp.dm8.act_sens = 3e-1; %--Change in oomph (E-field sensitivity) of DM8 actuators. Chosen empirically based on how much DM8 actuates during a control step. Larger than 1 is way too strong.

%%
mp.controller = 'plannedEFC';


% % PLANNED SEARCH EFC DEFAULTS     
mp.dm_ind = [1 2 8 9]; % vector of DMs used in controller at ANY time (not necessarily all at once or all the time). 
mp.ctrl.dmfacVec = 1;
%--CONTROL SCHEDULE. Columns of mp.ctrl.sched_mat are: 
    % Column 1: # of iterations, 
    % Column 2: log10(regularization), 
    % Column 3: which DMs to use (12, 128, 129, or 1289) for control
    % Column 4: flag (0 = false, 1 = true), whether to re-linearize
    %   at that iteration.
    % Column 5: flag (0 = false, 1 = true), whether to perform an
    %   EFC parameter grid search to find the set giving the best
    %   contrast .
    % The imaginary part of the log10(regularization) in column 2 is
    %  replaced for that iteration with the optimal log10(regularization)
    % A row starting with [0, 0, 0, 1...] is for relinearizing only at that time


SetJ = [...
    repmat([1,-5,129,1,1],[4,1]);...
    repmat([1,1j-1,129,1,1],[4,1]);...
    repmat([1,1j,129,1,1],[2,1]);...
    ];

mp.ctrl.sched_mat = [...
    repmat([1,1j,  128,1,1],[20,1]);...
    ...repmat([1,1j-1,12,1,1],[6,1]);...
    ...repmat([1,1j,  12,1,1],[1,1]);...
    repmat(SetJ,[10,1]);...
    ];

[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);


%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);

