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
close all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


%% Step 2: Load default model parameters

EXAMPLE_defaults_WFIRST_HLC_BMC_WFE

%% Properties for BMC Analysis

mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 6;           %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control

mp.d_dm1_dm2 = 0.6; % distance between DM1 and DM2 (meters)

%%--Pupil Masks
mp.P1.compact.Nbeam = 927;
mp.P2.compact.Nbeam = mp.P1.compact.Nbeam ;
mp.P3.compact.Nbeam = mp.P1.compact.Nbeam ;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam ;

mp.P1.full.Nbeam = 927; %--Number of pixel widths across the actual diameter of the beam/aperture (independent of beam centering)
mp.P2.full.Nbeam = mp.P1.full.Nbeam;
mp.P3.full.Nbeam = mp.P1.full.Nbeam;
mp.P4.full.Nbeam = mp.P1.full.Nbeam;

%--WFE maps and stops on DMs:
mp.dm1.inf_fn = 'influence_BMC_2kDM_400micron_res10.fits';
mp.dm2.inf_fn = 'influence_BMC_2kDM_400micron_res10.fits';

mp.flagDMwfe = true;
if(mp.flagDMwfe)
    mp.dm1.wfe = fitsread(sprintf('~/Data/BMC/wfe/bmc50_dm_wfe_%dpix_pupil.fits',mp.P1.full.Nbeam));
    mp.dm2.wfe = mp.dm1.wfe; %rot90(mp.dm1.wfe,2);
end

% mp.dm1.dm_spacing = 400e-6; %--User defined actuator pitch
% mp.dm2.dm_spacing = 400e-6; %--User defined actuator pitch

mp.flagDM1stop = false; %--logical flag whether to include the stop at DM1 or not
mp.flagDM2stop = true; %--logical flag whether to include the stop at DM2 or not
mp.dm2.Dstop = 49*mp.dm2.dm_spacing;  %--diameter of circular stop at DM2 and centered on the beam

mp.P2.D = 46.2937*mp.dm1.dm_spacing; % beam diameter at pupil closest to the DMs  (meters)
mp.P3.D = mp.P2.D;
mp.P4.D = mp.P2.D;

%mp.dm1.Nact = 50; % number of actuators across DM1
%mp.dm2.Nact = 50; % number of actuators across DM2

%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = false;%true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 36;
mp.TrialNum = 1;

%--Control
mp.Nitr = 50; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

%% [OPTIONAL] Start from a previous FALCO trial's DM settings

% fn_prev = sprintf('Series0033_Trial%04d_HLC_WFIRST180718_3DM50_z%s_IWA2.7_OWA10_6lams575nm_BW10_plannedEFC_snippet.mat',13,num2str(mp.d_dm1_dm2));
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);

