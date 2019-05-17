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
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path
% rmpath([mp.path.cvx 'lib/narginchk_:']) %--Legend plotting issue if CVX's narginchk function is used instead of Matlab's default function.

%% Step 2: Load default model parameters

EXAMPLE_defaults_LUVOIRB_VC_design

%% Step 3: Overwrite default values as desired

%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;
mp.flagFiber = true;

%--Record Keeping
mp.SeriesNum = 6482648;
mp.TrialNum = 1;

mp.fracBW = 0.5;
mp.Nsbp = 22;
mp.Nitr = 10;

%--LSB settings
mp.dm1.HminStep = 5e-11;
mp.dm2.HminStep = 5e-11;

%--Grid- or Line-Search Settings
% mp.ctrl.log10regVec = -5:1/2:0; %--log10 of the regularization exponents (often called Beta values)

%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

% %--DEBUGGING
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

%--Special settings for fibers

if(mp.flagFiber)
    mp.Fend.FOV = 1.5;
    mp.Fend.res = 15; %Has to be much higher than normal to avoid checkerboarding/ringing when going to F5.

    %--Fiber properties
    mp.fiber.a = 1.22; %Radius of the fiber core in lambda_0/D
    mp.fiber.NA = 2.5e-7; %Numerical aperture of the fiber in bizarro units
    
    %--Fiber tip plane properties (i.e., focal plane of lenslet(s)
    mp.F5.res = 4;
    mp.F5.FOV = 10;
    mp.F5.fiberPos = [0 0]; %Position of the fiber center in F5 in lambda/D.  
                            %Should be zero unless testing fiber/lenslet misalignments.
    
    %--Lenslet properties
    mp.Fend.lensletWavRad = 1.22; %Radius of the lenslet(s) in lambda_0/D
    mp.Fend.Nlens = 5; %Number of lenslets in Fend
    mp.Fend.x_lenslet = [4 9 14 19 24]; %Lenslet positions in Fend in lambda_0/D
    mp.Fend.y_lenslet = [0 0 0 0 0];
    mp.lensletFL = 150e-6; %Lenslet focal length in meters
    
    %--Off-axis, incoherent point source (exoplanet)
    mp.c_planet = 1e-10; %contrast of exoplanet
    mp.x_planet = -mp.Fend.x_lenslet(1); %Position of the exoplanet in lambda_0/D
    mp.y_planet = -mp.Fend.y_lenslet(1);
    %Note that the above coordinates are flipped from the lenslet positions
    %for some damn reason, so input the NEGATIVE of where you want the
    %planet to be.
    
    mp.thput_eval_x = mp.x_planet;
    mp.thput_eval_y = mp.y_planet;
end

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);