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

EXAMPLE_defaults_WFIRST_PhaseB_PROPER_HLC


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

%--GRID SEARCH EFC    
mp.controller = 'gridsearchEFC';
mp.Nitr = 5; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1;%1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.ctrl.flagUseModel = true; %--Whether to perform a model-based (vs empirical) grid search for the controller
mp.ctrl.log10regVec = -5:1:-2; %--log10 of the regularization exponents (often called Beta values)

%%

% mp.layout = 'wfirst_phaseb_simple';  %--Which optical layout to use. 'wfirst_phaseb_proper' or 'wfirst_phaseb_simple'


mp.Nsbp = 1;
mp.fracBW = 0.01;
% mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

% mp.full.pol_conds = 10; %--DEBUGGING ONLY
% mp.Nsbp = 3;
% mp.Nwpsbp = 3;%1;%3;%7;
% mp.fracBW = 0.10;
% mp.flagParfor = true; %--whether to use parfor for Jacobian calculation


% mp.Nsbp = 9;
% mp.fracBW = 0.10;%0.01;%0.04;
% mp.flagParfor = true; %--whether to use parfor for Jacobian calculation



% mp.Nsbp = 3;
% mp.fracBW = 0.025;
% mp.flagParfor = true; %--whether to use parfor for Jacobian calculation

% mp.Nsbp = 5;
% mp.fracBW = 0.05;
% mp.flagParfor = true; %--whether to use parfor for Jacobian calculation

if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

%%
lam_occ = lambdaFacs*mp.lambda0;

mp.F3.compact.Nxi = 40; mp.F3.compact.Neta = mp.F3.compact.Nxi;
mp.compact.FPMcube = zeros(mp.F3.compact.Nxi,mp.F3.compact.Nxi,mp.Nsbp);

prefix = '/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/hlc_20190210/run461_nro_';
fpm_axis = 'p';

for si=1:mp.Nsbp
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);
    
    %--Unclear which is the correct orientation yet
    fn_p_r = [prefix  'occ_lam' num2str(lam_occ(si),12) 'theta6.69pol'   fpm_axis   '_' 'real_crop.fits'];
    fn_p_i = [prefix  'occ_lam' num2str(lam_occ(si),12) 'theta6.69pol'   fpm_axis   '_' 'imag_crop.fits'];
%     fn_p_r = [prefix  'occ_lam' num2str(lam_occ(si),12) 'theta6.69pol'   fpm_axis   '_' 'real_rotated_crop.fits'];
%     fn_p_i = [prefix  'occ_lam' num2str(lam_occ(si),12) 'theta6.69pol'   fpm_axis   '_' 'imag_rotated_crop.fits'];
   
    mp.compact.FPMcube(:,:,si) = complex(fitsread(fn_p_r),fitsread(fn_p_i));


end
%%
for si=1:mp.Nsbp
   figure(100); imagesc(abs(mp.compact.FPMcube(:,:,si))); axis xy equal tight; colorbar; drawnow; pause(0.1); 
end
% return
%% Step 3b: Obtain the phase retrieval phase.

mp.full.input_field_rootname = '/Users/ajriggs/Repos/falco-matlab/data/maps/input_full';



optval.phaseb_dir = mp.full.phaseb_dir;

optval.cor_type = mp.full.cor_type;

optval.source_x_offset =0;
% optval.zindex = [];
% optval.zval_m = [];
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.use_errors = mp.full.use_errors;
optval.polaxis = mp.full.polaxis; 
% 
% % 2. full model, for regular psf 
% % EE = prop_run_multi(['wfirst_phaseb_v2'], lam_array, npsf, 'quiet', 'passvalue',optval );
% EE = prop_run(['wfirst_phaseb_v2'], lam_array, npsf, 'quiet', 'passvalue',optval );

% 3. full model, for field

optval.dm1_m = fitsread([mp.full.phaseb_dir 'dm1_flatten_pol10_575nm.fits']);
optval.use_dm1 =1;

optval.end_at_fpm_exit_pupil = 1;
optval.output_field_rootname = ['fld_at_xtPup'];
optval.use_fpm=0;
optval.use_hlc_dm_patterns=0;
nout = 1024;%512; 			% nout > pupil_daim_pix
% if testcase >2; nout =1024; end % >= pupil_daim_pix; 
% fld = prop_run_multi(['wfirst_phaseb_v2'], lam_array, nout, 'quiet', 'passvalue',optval );




mp.P1.compact.E = ones(ceil_even(mp.P1.compact.Nbeam+1),ceil_even(mp.P1.compact.Nbeam+1),mp.Nsbp); %--Initialize
for si=1:mp.Nsbp
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);

    fld = prop_run(['wfirst_phaseb_v2b'], lambda_um, nout, 'quiet', 'passvalue',optval );
    % % % fld(2:end,2:end) = rot90(fld(2:end,2:end),2);

    % figure(601); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv;
    % figure(602); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula;
    figure(605); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv;
    figure(606); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula;

    lams = num2str(lambda_um, '%6.4f');
    polaxis = 0;
    pols = ['polaxis'  num2str(polaxis,2)];
    fitswrite(real(fld), [mp.full.input_field_rootname '_' lams 'um_' pols '_real.fits' ]);
    fitswrite(imag(fld), [mp.full.input_field_rootname '_' lams 'um_' pols '_imag.fits' ]);

    %%--Downsampling for the compact model
    dxF = 1;
    dxC = mp.P1.full.Nbeam/mp.P1.compact.Nbeam;

    Nf = length(fld);
    Nc = ceil_even( (mp.P1.compact.Nbeam/mp.P1.full.Nbeam)*Nf );

    xF = (-Nf/2:Nf/2-1)*dxF;
    xC = (-Nc/2:Nc/2-1)*dxC;

    [Xf,Yf] = meshgrid(xF);
    [Xc,Yc] = meshgrid(xC);

    fldC = interp2(Xf,Yf,fld,Xc,Yc,'cubic',0); %--Downsample by interpolation

    figure(607); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv;
    figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula;

    %--Subtract out the tip/tilt
    Nbeam = 309;
    Narray = ceil_even(mp.P1.compact.Nbeam+1);

    x = (-Narray/2:Narray/2-1)/Nbeam; 
    [tip,tilt] = meshgrid(x);
    
    fldC = padOrCropEven(fldC,ceil_even(mp.P1.compact.Nbeam+1));
    fldC = fldC.*exp(1i*tilt*1.24*(1/lambdaFacs(si)));
    
    temp = 0*fldC;
    temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
    mp.P1.compact.E(:,:,si) = temp;
    
    figure(617+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv;


end

return
%% Find the coefficient of the tip to use on the phase retrieval map to flatten it
% Nbeam = 309;
% Narray = ceil_even(mp.P1.compact.Nbeam+1);
% 
% x = (-Narray/2:Narray/2-1)/Nbeam; 
% [tip,tilt] = meshgrid(x);
% % [tilt,tip] = meshgrid(x);
% 
% fldCcrop = padOrCropEven(fldC,Narray);
% ang = angle(fldCcrop);
% absF = abs(fldCcrop);
% mask = zeros(Narray);
% mask(absF>=0.5*max(absF(:))) = 1;
% mask_ele = find(mask==1);
% figure(610); imagesc(mask); axis xy equal tight; colorbar; colormap jet;
% figure(611); imagesc(mask.*ang); axis xy equal tight; colorbar; colormap jet;
% 
% 
% a2 = tip(:).'*ang(:);
% a3 = tilt(:).'*ang(:);
% piston = ones(size(tip));
% 
% facs = -1:0.01:1.5; %1.01:0.01:1.5;
% Nfacs = length(facs);
% rms_vec = zeros(Nfacs,1);
% 
% for ii=1:Nfacs
%     angNew = angle(fldCcrop.*exp(1i*tilt*facs(ii)));
%     rms_vec(ii) = falco_rms(angNew(mask_ele));
%     
% end
% figure(600); plot(facs,rms_vec);
% 
% % figure(610); imagesc(mask.*angle(fldCcrop.*exp(1i*tilt*1.2)) ); axis xy equal tight; colorbar; colormap jet;

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
