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

defaults_WFIRST_PhaseB_PROPER_SPC_Spec_custom_PIQ

%% Step 3a: Overwrite default values as desired

% %%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;
% mp.propMethodPTP = 'mft';

%--Record Keeping
mp.SeriesNum = 61;
mp.TrialNum = 2;%1;

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'Series...snippet.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

% %--DEBUGGING:
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
% mp.full.pol_conds = 10;

%% Step 3b: Obtain the phase retrieval phase.

mp.full.input_field_rootname = '/Users/ajriggs/Repos/falco-matlab/data/maps/input_full';
optval = mp.full;
optval.source_x_offset =0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.dm1_m = mp.full.dm1.flatmap;
optval.use_dm1 = 1;
optval.end_at_fpm_exit_pupil = 1;
optval.output_field_rootname = [fileparts(mp.full.input_field_rootname) filesep 'fld_at_xtPup'];
optval.use_fpm = 0;
optval.use_hlc_dm_patterns = 0;
nout = 1024; %512; 			% nout > pupil_daim_pix
optval.output_dim = 1024;%% Get the Input Pupil's E-field

optval.use_pupil_mask = false;  % No SP for getting initial phase

if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

mp.P1.compact.E = ones(mp.P1.compact.Nbeam+2,mp.P1.compact.Nbeam+2,mp.Nsbp); %--Initialize
for si=1:mp.Nsbp
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);

    fldFull = prop_run('model_full_wfirst_phaseb', lambda_um, nout, 'quiet', 'passvalue',optval );
    if(mp.flagPlot)
        figure(605); imagesc(angle(fldFull)); axis xy equal tight; colorbar; colormap hsv;
        figure(606); imagesc(abs(fldFull)); axis xy equal tight; colorbar; colormap parula;
    end

    lams = num2str(lambda_um, '%6.4f');
    polaxis = optval.polaxis;
    pols = ['polaxis'  num2str(polaxis,2)];
    fitswrite(real(fldFull), [mp.full.input_field_rootname '_' lams 'um_' pols '_real.fits' ]);
    fitswrite(imag(fldFull), [mp.full.input_field_rootname '_' lams 'um_' pols '_imag.fits' ]);

    %%--Downsampling for the compact model
    dxF = 1;
    dxC = mp.P1.full.Nbeam/mp.P1.compact.Nbeam;
    Nf = length(fldFull); %--N full
    Nc = ceil_even( (mp.P1.compact.Nbeam/mp.P1.full.Nbeam)*Nf ); %--N compact
    xF = (-Nf/2:Nf/2-1)*dxF;
    xC = (-Nc/2:Nc/2-1)*dxC;
    [Xf,Yf] = meshgrid(xF);
    [Xc,Yc] = meshgrid(xC);
    fldC = interp2(Xf,Yf,fldFull,Xc,Yc,'cubic',0); %--Downsample by interpolation
    fldC = padOrCropEven(fldC,ceil_even(mp.P1.compact.Nbeam+1));
    if(mp.flagPlot)
        figure(607); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula;
        figure(608+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv;
    end

    %--Assign to initial E-field in compact model.
    temp = 0*fldC;
    temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
    mp.P1.compact.E(:,:,si) = temp;
    
end

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);

