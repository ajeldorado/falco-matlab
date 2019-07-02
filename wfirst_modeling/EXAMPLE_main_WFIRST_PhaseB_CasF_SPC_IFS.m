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


%% Step 1: Define Necessary Paths on Your Computer System

%--Functions for when the full model uses PROPER
addpath('~/Repos/proper-models/wfirst_phaseb/matlab');

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

EXAMPLE_defaults_WFIRST_PhaseB_CasF_SPC_IFS


%% Step 3: Overwrite default values as desired


% %%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;
% mp.propMethodPTP = 'mft';

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

% %--DEBUGGING:
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

% mp.controller = 'plannedEFC';
% mp.ctrl.sched_mat = [...
%     repmat([1,1j,  12,1,1],[1,1]);...
%     repmat([1,1j-1,12,1,1],[10,1]);...
%     repmat([1,1j,  12,1,1],[1,1]);...
%     ];
% [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

% % % GRID SEARCH EFC    
mp.controller = 'gridsearchEFC';
mp.Nitr = 5; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

%% Step 3b: Obtain the phase retrieval phase.

mp.full.input_field_rootname = '/home/ajriggs/Repos/falco-matlab/data/maps/input_full';

optval.data_dir = mp.full.data_dir;

optval.cor_type = mp.full.cor_type;

optval.source_x_offset = 0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.use_errors = true;
optval.polaxis = 10; 
% 
% % 2. full model, for regular psf 
% % EE = prop_run_multi(['wfirst_phaseb_v2'], lam_array, npsf, 'quiet', 'passvalue',optval );
% EE = prop_run(['wfirst_phaseb_v2'], lam_array, npsf, 'quiet', 'passvalue',optval );

% 3. full model, for field

optval.dm1_m = fitsread ([mp.full.data_dir 'errors_polaxis10_dm.fits']);
optval.use_dm1 =1;

optval.end_at_fpm_exit_pupil =1;
optval.output_field_rootname = [fileparts(mp.full.input_field_rootname) filesep 'fld_at_xtPup'];
optval.use_fpm = 0;
optval.use_hlc_dm_patterns=0;
nout = 1024;%512; 			% nout > pupil_daim_pix
% if testcase >2; nout =1024; end % >= pupil_daim_pix; 
% fld = prop_run_multi(['wfirst_phaseb_v2'], lam_array, nout, 'quiet', 'passvalue',optval );

lambda_um = 1e6*mp.lambda0;

fld = prop_run(['model_full_wfirst_phaseb'], lambda_um, nout, 'quiet', 'passvalue',optval );
% % % fld(2:end,2:end) = rot90(fld(2:end,2:end),2);

% figure(601); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv;
% figure(602); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula;
figure(605); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv; drawnow;
figure(606); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula; drawnow;

lams = num2str(lambda_um, '%6.4f');
pols = ['polaxis'  num2str(optval.polaxis,2)];
fitswrite(real(fld), [mp.full.input_field_rootname '_' lams 'um_' pols '_real.fits' ]);
fitswrite(imag(fld), [mp.full.input_field_rootname '_' lams 'um_' pols '_imag.fits' ]);

%
dxF = 1;
dxC = mp.P1.full.Nbeam/mp.P1.compact.Nbeam;

Nf = length(fld);
Nc = ceil_even( (mp.P1.compact.Nbeam/mp.P1.full.Nbeam)*Nf );

xF = (-Nf/2:Nf/2-1)*dxF;
xC = (-Nc/2:Nc/2-1)*dxC;

[Xf,Yf] = meshgrid(xF);
[Xc,Yc] = meshgrid(xC);

fldC = interp2(Xf,Yf,fld,Xc,Yc,'cubic',0); %--Downsample by interpolation

figure(607); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula; drawnow;

% fldC = padOrCropEven(fldC,mp.P1.compact.Nbeam+2);
% temp = 0*fldC;
% temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
% mp.P1.compact.E = temp;

    fldC = padOrCropEven(fldC,ceil_even(mp.P1.compact.Nbeam+1));
    temp = 0*fldC;
    temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
    mp.P1.compact.E = temp; %mp.P1.compact.E(:,:,si) = temp;
    
    figure(617); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv;


%%
% return

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);
