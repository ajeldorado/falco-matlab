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

EXAMPLE_defaults_WFIRST_PhaseB_PROPER_SPC_IFS


%% Step 3: Overwrite default values as desired

% mp.full.output_dim = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
% mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D

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
% mp.fracBW = 0.08;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 3;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.flagParfor = true; %--whether to use parfor for Jacobian calculation

%--DEBUGGING:
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
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
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.ctrl.flagUseModel = true; %--Whether to perform a model-based (vs empirical) grid search for the controller


%% Step 3b: Obtain the phase retrieval phase.

mp.full.input_field_rootname = '/Users/ajriggs/Repos/falco-matlab/data/maps/input_full';


optval.phaseb_dir = mp.full.phaseb_dir;

optval.cor_type = mp.full.cor_type;

optval.source_x_offset =0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;

optval.use_errors = mp.full.use_errors;
optval.polaxis = mp.full.polaxis; 

optval.dm1_m = fitsread([mp.full.phaseb_dir 'dm1_flatten_pol10_730nm.fits']);
optval.use_dm1 = 1;

optval.end_at_fpm_exit_pupil = 1;
optval.output_field_rootname = ['fld_at_xtPup'];
optval.use_fpm=0;
optval.use_hlc_dm_patterns=0;
nout = 1024;%512; 			% nout > pupil_daim_pix
% if testcase >2; nout =1024; end % >= pupil_daim_pix; 
% fld = prop_run_multi(['wfirst_phaseb_v2'], lam_array, nout, 'quiet', 'passvalue',optval );


%% Get the Input Pupil's E-field

if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

mp.P1.compact.E = ones(mp.P1.compact.Nbeam+2,mp.P1.compact.Nbeam+2,mp.Nsbp); %--Initialize
for si=1:mp.Nsbp
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);

    fld = prop_run(['model_full_wfirst_phaseb'], lambda_um, nout, 'quiet', 'passvalue',optval );
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

%     %%--Downsampling for the compact model
%     dxF = 1;
%     dxC = mp.P1.full.Nbeam/mp.P1.compact.Nbeam;
% 
%     Nf = length(fld);
%     Nc = ceil_even( (mp.P1.compact.Nbeam/mp.P1.full.Nbeam)*Nf );
% 
%     xF = (-Nf/2:Nf/2-1)*dxF;
%     xC = (-Nc/2:Nc/2-1)*dxC;
% 
%     [Xf,Yf] = meshgrid(xF);
%     [Xc,Yc] = meshgrid(xC);
% 
%     fldC = interp2(Xf,Yf,fld,Xc,Yc,'cubic',0); %--Downsample by interpolation
% 
%     figure(607); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv;
%     figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula;
% 
% 
%     fldC = padOrCropEven(fldC,mp.P1.compact.Nbeam+2);
%     temp = 0*fldC;
%     temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
%     mp.P1.compact.E(:,:,si) = temp;
    
    
    
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

    %--Divide out the tip/tilt
    Narray = ceil_even(mp.P1.compact.Nbeam+1);

    x = (-Narray/2:Narray/2-1)/mp.P1.compact.Nbeam; 
    [tip,tilt] = meshgrid(x);
    
    fldC = padOrCropEven(fldC,ceil_even(mp.P1.compact.Nbeam+1));
    fldC = fldC.*exp(1i*tilt*0.87*(1/lambdaFacs(si)));
    
    temp = 0*fldC;
    temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
    mp.P1.compact.E(:,:,si) = temp;
    
    figure(617+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
    
    

end

%% Find the coefficient of the tip to use on the phase retrieval map to flatten it
% Narray = ceil_even(mp.P1.compact.Nbeam+1);
% 
% x = (-Narray/2:Narray/2-1)/mp.P1.compact.Nbeam; 
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


%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);
