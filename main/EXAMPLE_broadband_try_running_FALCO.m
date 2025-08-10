% Copyright 2018-2021 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to run wavefront control on a vortex coronagraph with a segmented
% input pupil.

clear
T_tot_start = tic;

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
mp.flagSaveWS = true;  %--Whether to save out entire (large) workspace at the end of trial. Default is false


%% Step 2: Load default model parameters

EXAMPLE_broadband_defaults_try_running_FALCO %*import and run the script


%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation %* false to disable parallel processing
mp.flagPlot = false;    %* enable plotting
mp.flag_timeMaya = false; 

%--Record Keeping
mp.SeriesNum = 11;
mp.TrialNum = 1;

%--Use just 1 wavelength for initial debugging/testing of code
mp.fracBW = 0.1;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 3;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
mp.Nsbp_bb = 3;    %mc number of sub-bandpasses for broadband pwp

mp.Nitr =10; %--Number of wavefront control iterations

mp.est.flagUseJac = true;
%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',...
    num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp); % Flesh out the rest of the workspace prior to wavefront estimation and control.
mkdir([mp.path.config filesep mp.runLabel]) 

[mp, out] = falco_wfsc_loop(mp, out);

T_tot = toc(T_tot_start)
