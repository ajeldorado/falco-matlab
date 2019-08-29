% Copyright 2019, B. Nemati and S. Miller. 
% -------------------------------------------------------------------------
% FRNcalculator Version_2.0 (updated 5/13/19)
% 
% This version is not for evaluating various detectors but for evaluating
% the yield with the nominal parameters
% IMG Mode


% clc; clear; close all;

% filepath = mfilename('fullpath');
% this_path = fileparts(filepath);
% parent_path = fileparts(this_path);
% addpath(fullfile(parent_path, 'calculateFRN')); % add full path of calculateFRN directory

addpath('~/Repos/falco-matlab/wfirst_modeling/FRN/calculateFRN/'); % add full path of calculateFRN directory

% Specify full path of folder contaning coronagraph data. See README to find what what
% files this directory requires.
CG_Directory = '~/Repos/falco-matlab/data/frn/CGdata_EB_IMG_NF_Band1';
% CG_Directory = fullfile(parent_path, 's_EB_IMG_NF_Band1/CGdata');
MUFcaseInput = 'Standard'; % choices are 'Standard' and 'Unity'
NItoContrast = 0.922; % user needs to calculate this for this annular zone
mode = 'IMG'; % choices are: 'IMG' 'IFS'

% scenarioName only matters if it is one of the following:
%   'EB IMG NF - Band 1',
%   'EB SPEC IFS - Band 3
%   'EB IMG WF = Band 4'
scenarioName  = 'EB IMG NF - Band 1';
detType       = 'EM_PC_SRR';
CGtauPol      = 1; % (1 or 0.5): if 0.5, the design file includes polarizer loss
centerLambda  = 575 * 1e-9; % meters
bandWidth     = 0.10;
yearsL2       = 63/12; % for detector
SNR           = 10; % SNR desired
frameTime     = 6; % seconds
k_pp          = 2;
maxIntegTime  = 10 * 3600; % seconds
dutyFactor    = 0.70;

% set to true to suppress text output from calculateFRN
suppressOutput = false;

% Choose planet from planet table. To see planet values use:
% readPlanetTable(planetNo)
% Band numbers are only supported for planet number 1 (EB Fiducial Planet).
% Supported band numbers are 1, 3, and 4.
% Band number other than 1 will alter planet.A and planet.Rp_R_j
planetNo = 1;
bandNo = 1;

tic
FRN = calculateFRN(CG_Directory, MUFcaseInput, NItoContrast, mode, scenarioName, detType, CGtauPol, centerLambda, bandWidth, yearsL2, SNR, frameTime, k_pp, maxIntegTime, dutyFactor, suppressOutput, planetNo, bandNo);

fprintf('\n\nFRN Values\n');
fprintf('--------------------------\n');
fprintf('Photometry         = %.3f\n', FRN.photometry);
fprintf('Contrast Stability = %.3f\n', FRN.cstability);
fprintf('Total              = %.3f\n\n', FRN.total);
toc
