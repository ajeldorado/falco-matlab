%%% Paramameters Function
%
% This function will be called by test scripts requiring one or more of the
% mp parameters necessary to run the test for given FALCO function/script.
% No input is required.

function [mp] = Parameters()
%% Define Necessary Paths on Your Computer System

% In this section we define and add necessary paths to FALCO.
mp.path.falco = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))); % falco-matlab directory;
addpath(genpath([mp.path.falco filesep 'setup']))
addpath(genpath([mp.path.falco filesep 'lib']))
addpath(genpath([mp.path.falco filesep 'config']))

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

%%% Use just 1 wavelength for initial debugging of code
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1; 

%% Generate the label associated with this trial
%
% This generates the label associated with a particular trial, not really
% needed here.
mp.runLabel = 'testing_label';


%% Perform the Wavefront Sensing and Control
%
% This section generates parameters required to run the full wfsc
% calculation.
[mp, out] = falco_flesh_out_workspace(mp);
% timerValEnd = toc;
% disp(['Time Elapsed:  ' num2str(timerValEnd)])
end
