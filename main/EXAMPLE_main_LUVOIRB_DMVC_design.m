% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform a DMVC simple design run.
%  1) Load the default model parameters for a vortex.
%  2) Specify the values to overwrite.
%  3) Run a single trial of WFC using FALCO.
%
% REVISION HISTORY:
% --------------
% Modified on 2019-02-26 by A.J. Riggs to load the defaults first.
% ---------------

clear all;
close all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = 'C:\CoronagraphSims\falco-matlab\';  %--Location of FALCO
mp.path.proper = 'C:\CoronagraphSims\FALCO\lib\PROPER\'; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = 'C:\CoronagraphSims\falco-matlab\data\brief\'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = 'C:\CoronagraphSims\falco-matlab\data\ws\'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


%% Step 2: Load default model parameters

EXAMPLE_defaults_LUVOIRB_VC_design


%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%%-- segmented mirror errors
mp.numSegments = hexSegMirror_numSegments(4); % Number of segments in "full" hex aperture
% LUVOIR B has four rings, but ignores some corner segmentes 

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
fn_prev = 'Series0867_Trial5309_Vortex_LUVOIR_B_offaxis_2DM64_z0.8_IWA2_OWA26_1lams400nm_BW2.5_gridsearchEFC_snippet.mat';
temp = load(fn_prev,'out');
mp.dm1.V = temp.out.DM1V;
mp.dm2.V = temp.out.DM2V;
clear temp

% %--DEBUGGING
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.lambda0 = 400e-9;
% mp.flagParfor = false; %--whether to use parfor for Jacobian calculation


%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);

