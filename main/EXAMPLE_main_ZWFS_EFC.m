% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

clear;

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations
mp.path.falco = 'C:\CoronagraphSims\falco-matlab\';  %--Location of FALCO
mp.path.proper = 'C:\CoronagraphSims\FALCO\lib\PROPER\'; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = [mp.path.falco,'data/brief/']; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = [mp.path.falco,'data/ws/']; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%---Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

%% Step 2: Load default model parameters

EXAMPLE_defaults_LUVOIRB_VC_design

%% Step 3: Overwrite default values as desired

mp.SeriesNum = 1;
mp.TrialNum = 99;

mp.flagParfor = false;
mp.flagPlot = true;
mp.flagFiber = false;
mp.flagLenslet = false;
mp.flagZWFS = true;

%--[OPTIONAL] Start from a previous FALCO trial's DM settings
fn_prev = 'Series0867_Trial5309_Vortex_LUVOIR_B_offaxis_2DM64_z0.8_IWA2_OWA26_1lams400nm_BW2.5_gridsearchEFC_snippet.mat';
temp = load(fn_prev,'out');
mp.dm1.V = temp.out.DM1V;
mp.dm2.V = temp.out.DM2V;
clear temp

%%--Coronagraph and Pupil Type
mp.coro = 'Vortex';    %--Tested Options: 'LC','HLC','SPLC','Vortex'
mp.flagApod = true;
mp.whichPupil = 'LUVOIR_B_offaxis';
mp.F3.VortexCharge = 6;

%%--Bandwidth and Wavelength Specs
mp.lambda0 = 400e-9; % central wavelength of bandpass (meters)
mp.fracBW = 10e-9/mp.lambda0;%0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 1;  % number of wavelengths or sub-bandpasses (sbp) across entire spectral band

%%-- segmented mirror errors
numSegments = hexSegMirror_numSegments(4); % Number of segments in "full" hex aperture
% LUVOIR B has four rings, but ignores some corner segmentes 

%%-- Focal Plane Mask (F3) Properties

pupilDiam_m = mp.P2.D;%m
maskDepth_m  = 213e-9;%m
maskMaterial = 'FS';%Fused Silica 
maskRadius_lamOverD = 1.06/2;%(maskDiam_m/2)/fNum/mp.lambda0;

mp.F3.Rin = maskRadius_lamOverD;
mp.F3.Rout = Inf;
mp.F3.t = maskDepth_m; % Depth of the Roddier/Zernike spot
mp.FPMampFac = 1.0;% Transmission of the Roddier/Zernike spot
mp.FPMmaterial = maskMaterial;

mp.F3.full.res = 100;
mp.F3.compact.res = 100;

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -8:1/2:-2; %--log10 of the regularization exponents (often called Beta values)

%% Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Perform wavefront sensing and control

out = falco_wfsc_loop(mp);