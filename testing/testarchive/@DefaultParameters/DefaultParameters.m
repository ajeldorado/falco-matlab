%%% Default Paramameters Function
%
% This function will be called by test scripts requiring the smallest subset
% of the FALCO parameters (default parameters). This is necessary to test
% the subfunctions called inside the function falco_flesh_out_workspace.m. 
% For example, when we want to test the subfunction
% falco_set_optional_variables.m we do not need to run the full function to
% generate the FALCO's parameter structure (Parameters.m generates the full 
% set of FALCO's parameters).

function [mp] = DefaultParameters()
%% Define Necessary Paths on Your Computer System
%
% In this section we define and add necessary paths to FAlCO and PROPER. If 
% we do not define and add these paths we will not be able to call FALCO or PROPER
% functions.
mp.path.falco = '../../../falco-matlab'; 

addpath(genpath(mp.path.falco)) 
mp.path.proper = '/Users/lmarchen/Documents/PROPER';
addpath(genpath(mp.path.proper)) 

%% Load default model parameters
%
% This section runs a script found in falco repo to generate certain FALCO
% default parameters. These parameters can be overwritten below.
EXAMPLE_defaults_WFIRST_LC

%% Overwrite default values as desired
%
% This section overwrites some of the default parameters needed to run
% FALCO's main example.

%%%--command line message
mp.flagfprintf = false;

%%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = false;
mp.propMethodPTP = 'mft';

%%% Record Keeping
mp.TrialNum = 1;
mp.SeriesNum = 1;

%%% Use just 1 wavelength for initial debugging of code
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control

mp.F3.Rin = 2.7;    % maximum radius of inner part of the focal plane mask [lambda0/D]
mp.F3.RinA = mp.F3.Rin;   % inner hard-edge radius of the focal plane mask [lambda0/D]. Needs to be <= mp.F3.Rin
mp.Fend.corr.Rin = mp.F3.Rin;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.score.Rin = mp.F3.Rin;  % inner radius of dark hole scoring region [lambda0/D]

mp.P4.IDnorm = 0.45; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.78; %--Lyot stop OD [Dtelescope]

%% Generate the label associated with this trial
%
% This generates the label associated with a particular trial, not really
% needed here.
mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

return
