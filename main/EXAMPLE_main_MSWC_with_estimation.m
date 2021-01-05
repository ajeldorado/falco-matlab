% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform a DMVC simple design run.
%  1) Load the default model parameters for a vortex.
%  2) Specify the values to overwrite.
%  3) Run a single trial of WFC using FALCO.


clear
% close all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/proper_v3.2_matlab_11feb20/'; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


%% Step 2: Load default model parameters

EXAMPLE_defaults_MSWC


%% Step 3: Overwrite default values as desired

mp.jac.star.weights = [1, 1]; % star weights for control

% % Both Stars
starWeights = [1, 1];
mp.compact.star.count = 2;
mp.compact.star.xiOffsetVec = [0, 50]; %[0, 56];
mp.compact.star.etaOffsetVec = [0, -6];
mp.compact.star.weights = starWeights; % relative stellar peak intensities

mp.star.count = 2;
mp.star.xiOffsetVec = [0, 50]; %[0, 56];
mp.star.etaOffsetVec = [0, -6];
mp.star.weights = starWeights; % relative stellar peak intensities
mp.Fend.xiFOV = 66;
mp.Fend.etaFOV = 10;


% % % SNWC Only
% starWeights = [1];
% mp.compact.star.count = 1;
% mp.compact.star.xiOffsetVec = [50]; %[0, 56];
% mp.compact.star.etaOffsetVec = [-6];
% mp.compact.star.weights = starWeights; % relative stellar peak intensities
% 
% mp.star.count = 1;
% mp.star.xiOffsetVec = [50]; %[0, 56];
% mp.star.etaOffsetVec = [-6];
% mp.star.weights = starWeights; % relative stellar peak intensities
% mp.Fend.xiFOV = 66;
% mp.Fend.etaFOV = 10;
% 
% %--New variables for pairwise probing estimation:
% mp.est.probe.Npairs = 2;     % Number of pair-wise probe PAIRS to use.
% mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
% mp.probe.est.width = [6]; % Width of probed rectangular region. Units of lambda/D.
% mp.probe.est.height = [5]; % Height of probed rectangular region. Units of lambda/D.
% mp.est.probe.xiOffset = [6]; % xi (horizontal) offset of probed region's center in focal plane. Units of lambda/D.
% mp.est.probe.etaOffset = [6]; % eta (horizontal) offset of probed region's center in focal plane. Units of lambda/D.
% mp.est.probe.xOffset = 0;   % offset of probe center in x at DM [actuators]. Use to avoid central obscurations.
% mp.est.probe.yOffset = 10;    % offset of probe center in y at DM [actuators]. Use to avoid central obscurations.
% mp.est.probe.gainFudge = [5];     % empirical fudge factor to make average probe amplitude match desired value.







%--Correction and scoring region definition
mp.Fend.corr.Rin = 0;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 2;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 0;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 2;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 180;  % angular opening of dark hole scoring region [degrees]

mp.Fend.shape = 'square';
mp.Fend.xiOffset = 6;

mp.ctrl.log10regVec = -6:1/2:0; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.flagUseModel = false;%true;

% mp.Fend.xiFOV = 40;
% mp.Fend.etaFOV = 12;


% mp.jac.zerns = [1,2,3];
% mp.jac.Zcoef = 1e-9*[1,1,1];

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 2;

%--Use just 1 wavelength for initial debugging of code
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
mp.fracBW = 0.03;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 3;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 3;          %--Number of wavelengths to used to approximate an image in each sub-bandpass



mp.Nitr = 20; %--Number of wavefront control iterations

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control


[mp, out] = falco_flesh_out_workspace(mp);

% Apply a grid of spots to the input pupil to allow SNWC
block0 = ones(5, 5);
block0(3, 3) = 0.7;%0;
dotGrid = repmat(block0, [51, 51]);
dotGrid = pad_crop(dotGrid, mp.P1.compact.Narr, 'extrapval', 1);
% figure; imagesc(dotGrid); axis xy equal tight; colorbar;
for si = 1:mp.Nsbp
    wvl = mp.sbp_centers(si);
    mp.P1.compact.E(:, :, si) = dotGrid .* ones(mp.P1.compact.Narr);
    for wi = 1:mp.Nwpsbp
        mp.P1.full.E(:, :, wi, si) = dotGrid .* ones(mp.P1.full.Narr);
    end
end

[mp, out] = falco_wfsc_loop(mp, out);

