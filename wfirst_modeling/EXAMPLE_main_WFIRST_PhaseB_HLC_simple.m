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

%--Functions and some FITS files for when the full model uses PROPER
addpath('~/Repos/proper-models/wfirst_cgi/models_phaseb/matlab');
addpath('~/Repos/proper-models/wfirst_cgi/models_phaseb/matlab/examples');

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

EXAMPLE_defaults_WFIRST_PhaseB_HLC_simple


%% Step 3: Overwrite default values as desired

% %%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;
% mp.propMethodPTP = 'mft';

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%--Apply the HLC DM shapes from the design
mp.dm1.V = fitsread('hlc_dm1.fits')./mp.dm1.VtoH;
mp.dm2.V = fitsread('hlc_dm2.fits')./mp.dm2.VtoH;

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

% %--DEBUGGING:
% mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 1; % 1, 3, or 7  %--Number of wavelengths to used to approximate an image in each sub-bandpass
% % mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

% %--PLANNED EFC
% mp.controller = 'plannedEFC';
% mp.ctrl.sched_mat = [...
%     repmat([1,1j,12,1,1],[5,1]);...
%     repmat([1,1j-1,12,1,1],[25,1]);...
%     repmat([1,1j,12,1,1],[1,1]);...
%     ];
% [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%--GRID SEARCH EFC    
mp.controller = 'gridsearchEFC';
mp.Nitr = 10; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.ctrl.flagUseModel = false; %--Whether to perform a model-based (vs empirical) grid search for the controller
mp.ctrl.log10regVec = -6:1/2:-1; %--log10 of the regularization exponents (often called Beta values)

%% FPM representation vs wavelength: compact model

if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

lam_occ = lambdaFacs*mp.lambda0;

mp.F3.compact.Nxi = 40; %--Crop down to minimum size of the spot
mp.F3.compact.Neta = mp.F3.compact.Nxi;
mp.compact.FPMcube = zeros(mp.F3.compact.Nxi,mp.F3.compact.Nxi,mp.Nsbp);
fpm_axis = 'p';

for si=1:mp.Nsbp
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);
    fn_p_r = [mp.full.data_dir 'hlc_20190210/run461_occ_lam' num2str(lam_occ(si),12) 'theta6.69pol'   fpm_axis   '_' 'real.fits'];
    fn_p_i = [mp.full.data_dir 'hlc_20190210/run461_occ_lam' num2str(lam_occ(si),12) 'theta6.69pol'   fpm_axis   '_' 'imag.fits'];   
    mp.compact.FPMcube(:,:,si) = padOrCropEven(complex(fitsread(fn_p_r),fitsread(fn_p_i)),mp.F3.compact.Nxi);
end

%%
for si=1:mp.Nsbp
   if(mp.flagPlot);  figure(91); imagesc(abs(mp.compact.FPMcube(:,:,si))); axis xy equal tight; colorbar; title(sprintf('%d',si)); drawnow; pause(0.2);  end
end

%% FPM representation vs wavelength: full model

if(mp.Nsbp==1 && mp.Nwpsbp==1)
    lambdaFacs = 1;
elseif(mp.Nsbp==1 && mp.Nwpsbp>1)
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nwpsbp);
elseif(mp.Nsbp>1 && mp.Nwpsbp==1)
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
elseif(mp.Nsbp>1 && mp.Nwpsbp>1)
    Nlam = mp.Nsbp*mp.Nwpsbp - (mp.Nsbp-1);
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,Nlam);    
end
Nlam = length(lambdaFacs);

lam_occ = lambdaFacs*mp.lambda0;
mp.F3.full.Nxi = mp.F3.compact.Nxi; %--Same size as for compact model
mp.F3.full.Neta = mp.F3.full.Nxi;

for ilam=1:Nlam
    lambda_um = 1e6*mp.lambda0*lambdaFacs(ilam);
    fn_p_r = [mp.full.data_dir 'hlc_20190210/run461_occ_lam' num2str(lam_occ(ilam),12) 'theta6.69pol'   fpm_axis   '_' 'real.fits'];
    fn_p_i = [mp.full.data_dir 'hlc_20190210/run461_occ_lam' num2str(lam_occ(ilam),12) 'theta6.69pol'   fpm_axis   '_' 'imag.fits'];
    mp.full.FPMcube(:,:,ilam) = padOrCropEven(complex(fitsread(fn_p_r),fitsread(fn_p_i)),mp.F3.full.Nxi); 
end

for ilam=1:Nlam
   if(mp.flagPlot);  figure(92); imagesc(abs(mp.full.FPMcube(:,:,ilam))); axis xy equal tight; colorbar; title(sprintf('%d',ilam)); drawnow; pause(0.2);  end
end


%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

[out,mp] = falco_wfsc_loop(mp);

