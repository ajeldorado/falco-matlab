% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to run high-order WFSC with any of the Roman CGI's mask configurations.
%
% Requires the Roman CGI Phase C model and data from https://sourceforge.net/projects/cgisim/
%
% Refer to Figure 2 in the paper at https://arxiv.org/abs/2108.05986 to see
% all the flight mask configurations that can be modeled with this script.

clear

%% Define Necessary Paths on Your Computer System

% %%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
% mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


%% Uncomment the config file for the mask configuration that you want

% %--Officially SUPPORTED mask configurations:
EXAMPLE_config_Roman_CGI_HLC_NFOV_Band1()
% EXAMPLE_config_Roman_CGI_SPC_Bowtie_Band2()
% EXAMPLE_config_Roman_CGI_SPC_Bowtie_Band3()
% EXAMPLE_config_Roman_CGI_SPC_WFOV_Band4()

% %--UNSUPPORTED but included mask configurations:
% EXAMPLE_config_Roman_CGI_SPC_RotatedBowtie_Band2()
% EXAMPLE_config_Roman_CGI_SPC_RotatedBowtie_Band3()
% EXAMPLE_config_Roman_CGI_HLC_NFOV_Band2()
% EXAMPLE_config_Roman_CGI_HLC_NFOV_Band3()
% EXAMPLE_config_Roman_CGI_HLC_NFOV_Band4()
% EXAMPLE_config_Roman_CGI_SPC_WFOV_Band1()
% EXAMPLE_config_Roman_CGI_SPC_Multistar_Band1()
% EXAMPLE_config_Roman_CGI_SPC_Multistar_Band4()

%% Overwrite default values as desired

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


%% SETTINGS FOR QUICK RUN: SINGLE WAVELENGTH, SINGLE POLARIZATION, AND NO PROBING

mp.fracBW = 0.01; %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1; %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1; %--Number of wavelengths to used to approximate an image in each sub-bandpass
mp.full.pol_conds = 10;% [-2,-1,1,2]; %--Which polarization states to use when creating an image.
mp.estimator = 'perfect';
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation


%% Keep only the central bandpass's FPM if using just one wavelength with HLC

if (mp.Nsbp == 1) && strcmpi(mp.coro, 'HLC')
    nSlices = size(mp.compact.FPMcube, 3);
    mp.compact.FPMcube = mp.compact.FPMcube(:, :, ceil(nSlices/2));
end


%% Perform an idealized phase retrieval (get the E-field directly) of the entrance pupil

optval = mp.full;
optval.source_x_offset = 0;
optval.use_dm1 = true;
optval.use_dm2 = true;
nout = 1024;
optval.output_dim = 1024;
optval.use_pupil_mask = false;  % No SPM for getting entrance pupil phase
optval.use_fpm = false;
optval.use_lyot_stop = false;
optval.use_field_stop = false;
optval.use_pupil_lens = true;
optval = rmfield(optval, 'final_sampling_lam0');

% Use non-SPC flat maps for SPC since SPM has separate aberrations
% downstream that can't be fully captured at entrance pupil with the SPM in
% place. The SPM aberrations are flattened in a separate step not included
% here.
is_spc = strfind(lower(mp.coro), 'sp');
if is_spc
    optval.dm1_m = mp.full.dm1.flatmapNoSPM;
    optval.dm2_m = mp.full.dm2.flatmapNoSPM;
else
    optval.dm1_m = mp.full.dm1.flatmap;
    optval.dm2_m = mp.full.dm2.flatmap;
end

if mp.Nsbp == 1
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2, 1+mp.fracBW/2, mp.Nsbp);
end

%--Get the Input Pupil's E-field and downsample for the compact model
nCompact = ceil_even(mp.P1.compact.Nbeam+1);
mp.P1.compact.E = ones(nCompact, nCompact, mp.Nsbp);
for iSubband = 1:mp.Nsbp
    
    lambda_um = 1e6*mp.lambda0*lambdaFacs(iSubband);
    
    % Get aberrations for the full optical train
    optval.pinhole_diam_m = 0; % 0 means don't use the pinhole at FPAM
    fieldFullAll = prop_run('roman_phasec', lambda_um, nout, 'quiet', 'passvalue', optval);
    
    % Put pinhole at FPM to get back-end optical aberrations
    optval.pinhole_diam_m = mp.F3.pinhole_diam_m;
    fieldFullBackEnd = prop_run('roman_phasec', lambda_um, nout, 'quiet', 'passvalue', optval);
    optval.pinhole_diam_m = 0; % 0 means don't use the pinhole at FPAM
    
    % Subtract off back-end phase aberrations from the phase retrieval estimate
    phFrontEnd = angle(fieldFullAll) - angle(fieldFullBackEnd);
    swMask = logical(ampthresh(fieldFullAll));
    [phFrontEnd, ~] = removeZernikes(phFrontEnd, [0 1 1], [0 1 -1], swMask); % Remove tip/tilt/piston
    
    % Put front-end E-field into compact model
    fieldFull = abs(fieldFullAll).*exp(1j*phFrontEnd);
    fieldCompact = falco_filtered_downsample(fieldFull, mp.P1.compact.Nbeam/mp.P1.full.Nbeam, mp.centering);
    fieldCompact = pad_crop(fieldCompact, nCompact);
    mp.P1.compact.E(:, :, iSubband) = propcustom_relay(fieldCompact, 1, mp.centering); % Assign to initial E-field in compact model

    if mp.flagPlot
        figure(607); imagesc(angle(fieldCompact)); axis xy equal tight; colorbar; colormap hsv; drawnow;
        figure(608); imagesc(abs(fieldCompact)); axis xy equal tight; colorbar; colormap parula; drawnow;
    end
    
end

% Don't double count the pupil amplitude with the phase retrieval and a model-based mask
mp.P1.compact.mask = ones(size(mp.P1.compact.mask));

%% Generate the label associated with this trial

mp.runLabel = sprintf('Series%04d_Trial%04d_Roman_CGI_SPC_WFOV_', mp.SeriesNum, mp.TrialNum);


%% Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);
