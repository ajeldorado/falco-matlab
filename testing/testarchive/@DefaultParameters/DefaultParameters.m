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
mp.runLabel = 'testing_label';

end
