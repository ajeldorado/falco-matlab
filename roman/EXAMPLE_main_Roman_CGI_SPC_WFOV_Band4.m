% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to run the high-order wavefront correction with the Roman CGI's
% SPC-WFOV mode for Phase C.
%
% Requires the Roman CGI Phase C model and data from https://sourceforge.net/projects/cgisim/

clear

%% Step 1: Define Necessary Paths on Your Computer System

%--Tell Matlab where to find the PROPER model prescription and FITS files
addpath(genpath('~/Documents/Sim/cgi/public/roman_phasec_v1.2.4/matlab/'));

% %%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
% mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


%% Step 2: Load default model parameters

EXAMPLE_config_Roman_CGI_SPC_WFOV_Band4()


%% Step 3: Overwrite default values as desired

% %%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'SeriesM_TrialN_snippet.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.dm1.Vall(:, :, end);
% mp.dm2.V = temp.out.dm2.Vall(:, :, end);
% clear temp

mp.controller = 'plannedEFC';
mp.ctrl.sched_mat = repmat([1, 1j, 12, 0, 1], [5, 1]);
% mp.ctrl.sched_mat = [...
%     [0,0,0,1,0];
%     repmat([1,1j,12,0,1],[5,1]);...
%     [1,-5,12,0,0];...
%     repmat([1,1j,12,0,1],[9,1]);...
%     ];
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%--QUICK CHECK WITH MONOCHROMATIC LIGHT
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

%% Step 3b: Perform an idealized phase retrieval (get the E-field directly)

optval = mp.full;
optval.source_x_offset = 0;
optval.use_dm1 = true;
optval.dm1_m = mp.full.dm1.flatmap;
optval.use_dm2 = true;
optval.dm2_m = mp.full.dm2.flatmap;
optval.end_at_fpm_exit_pupil = 1;
optval.use_fpm = 0;
nout = 1024;
optval.output_dim = 1024;
optval.use_pupil_mask = false;  % No SPM for getting initial phase

if mp.Nsbp == 1
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2, 1+mp.fracBW/2, mp.Nsbp);
end

%--Get the Input Pupil's E-field and downsample for the compact model
nCompact = ceil_even(mp.P1.compact.Nbeam+1);
mp.P1.compact.E = ones(nCompact, nCompact, mp.Nsbp);
for si = 1:mp.Nsbp
    
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);
    fldFull = prop_run('roman_phasec', lambda_um, nout, 'quiet', 'passvalue', optval);
    fldC = falco_filtered_downsample(fldFull, mp.P1.compact.Nbeam/mp.P1.full.Nbeam, mp.centering);
    fldC = pad_crop(fldC, nCompact);
    mp.P1.compact.E(:, :, si) = propcustom_relay(fldC, 1, mp.centering); % Assign to initial E-field in compact model

    if mp.flagPlot
        figure(607); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
        figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula; drawnow;
    end
    
end

% Don't double count the pupil amplitude with the phase retrieval and a model-based mask
mp.P1.compact.mask = ones(size(mp.P1.compact.mask));

%% Step 4: Generate the label associated with this trial

mp.runLabel = sprintf('Series%04d_Trial%04d_Roman_CGI_SPC_WFOV_', mp.SeriesNum, mp.TrialNum);


%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);
