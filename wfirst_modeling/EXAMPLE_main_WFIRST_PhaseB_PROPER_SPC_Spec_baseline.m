% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to run the wavefront correction with the Roman (formerly WFIRST)
% CGI's SPC-Spec mode from Phase B. 

clear


%% Step 1: Define Necessary Paths on Your Computer System

%--Tell Matlab where to find the PROPER model prescription and FITS files
addpath('~/Repos/proper-models/wfirst_cgi/models_phaseb/matlab');
addpath('~/Repos/proper-models/wfirst_cgi/models_phaseb/matlab/examples');

% %%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
% mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


%% Step 2: Load default model parameters

EXAMPLE_defaults_WFIRST_PhaseB_PROPER_SPC_Spec


%% Step 3: Overwrite default values as desired

% %%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 49;
mp.TrialNum = 2;

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

mp.controller = 'plannedEFC';
mp.ctrl.sched_mat = repmat([1,1j,12,0,1],[5,1]);
% mp.ctrl.sched_mat = [...
%     [0,0,0,1,0];
%     repmat([1,1j,12,0,1],[5,1]);...
%     [1,-5,12,0,0];...
%     repmat([1,1j,12,0,1],[9,1]);...
%     ];
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);


%% Choose starting DM commands to flatten the wavefront

mp.full.input_field_rootname = '/Users/ajriggs/Repos/falco-matlab/data/maps/input_full';
optval = mp.full;
optval.source_x_offset =0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.use_dm1 = false;
optval.use_dm2 = false;

optval.end_at_fpm_exit_pupil = 1;
optval.output_field_rootname = [fileparts(mp.full.input_field_rootname) filesep 'fld_at_xtPup'];
optval.use_fpm = 0;
optval.use_hlc_dm_patterns = 0;
nout = 1024; %512; 			% nout > pupil_daim_pix
optval.output_dim = 1024;%% Get the Input Pupil's E-field

optval.use_pupil_mask = false;  % No SP for getting initial phase

if mp.Nsbp == 1
    lambdaFacs = 1;
elseif mp.Nwpsbp == 1
    lambdaFacs = linspace(1-mp.fracBW/2, 1+mp.fracBW/2, mp.Nsbp);
else
    DeltaBW = mp.fracBW/(mp.Nsbp)*(mp.Nsbp-1)/2;
    lambdaFacs = linspace(1-DeltaBW, 1+DeltaBW, mp.Nsbp);
end

%--Get the Input Pupil's E-field
mp.P1.compact.E = ones(mp.P1.compact.Nbeam+2,mp.P1.compact.Nbeam+2,mp.Nsbp); %--Initialize
si = ceil(mp.Nsbp/2);
lambda_um = 1e6*mp.lambda0*lambdaFacs(si);

EpupFull = prop_run('model_full_wfirst_phaseb', lambda_um, nout, 'quiet', 'passvalue', optval );
if(mp.flagPlot)
    figure(505); imagesc(angle(EpupFull)); axis xy equal tight; colorbar; colormap hsv; drawnow;
    figure(506); imagesc(abs(EpupFull)); axis xy equal tight; colorbar; colormap parula; drawnow;
end

% Which pixels to use when flattening the phase
mask = 0*EpupFull;
mask(abs(EpupFull) > 1e-1*max(abs(EpupFull(:)))) = 1;
figure(3); imagesc(mask); axis xy equal tight; colorbar; drawnow;

    
surfaceToFit = -0.5*mask.*angle(EpupFull) * (mp.lambda0/(2*pi));
figure(1); imagesc(abs(EpupFull)); axis xy equal tight; colorbar; colormap parula; drawnow;
figure(2); imagesc(surfaceToFit); axis xy equal tight; colorbar; colormap parula;  drawnow;

mp.dm1.inf0 = fitsread(mp.dm1.inf_fn);
% mp.dm1.dm_spacing = 400e-6;
mp.dm1.dx_inf0 = mp.dm1.dm_spacing/10;
mp.dm1.dx = mp.P2.D/mp.P2.compact.Nbeam;
mp.dm1.centering = 'pixel';
dm1copy = mp.dm1;
dm1copy.dx = mp.dm1.dx * (mp.P1.compact.Nbeam/mp.P1.full.Nbeam);
flatmap = falco_fit_dm_surf(dm1copy, surfaceToFit);
figure(4); imagesc(flatmap); axis xy equal tight; colorbar; drawnow;

% Split commands evenly between DMs 1 and 2
mp.full.dm1.flatmap = flatmap/2;% ./ mp.dm1.VtoH;
mp.full.dm2.flatmap = flatmap/2;% ./ mp.dm2.VtoH;


%% Step 3b: Perform an idealized phase retrieval (get the E-field directly)

mp.full.input_field_rootname = '/Users/ajriggs/Repos/falco-matlab/data/maps/input_full';
optval = mp.full;
optval.source_x_offset =0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.use_dm1 = true;
optval.dm1_m = mp.full.dm1.flatmap;% .* mp.dm1.VtoH;
optval.use_dm2 = true;
optval.dm2_m = mp.full.dm2.flatmap;% .* mp.dm2.VtoH;

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

%--Get the Input Pupil's E-field
nCompact = ceil_even(mp.P1.compact.Nbeam+1);
mp.P1.compact.E = ones(nCompact, nCompact, mp.Nsbp);
for si = 1:mp.Nsbp
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);

    fldFull = prop_run('model_full_wfirst_phaseb', lambda_um, nout, 'quiet', 'passvalue',optval );
    if(mp.flagPlot)
        figure(605); imagesc(angle(fldFull)); axis xy equal tight; colorbar; colormap hsv; drawnow;
        figure(606); imagesc(abs(fldFull)); axis xy equal tight; colorbar; colormap parula; drawnow;
    end

    lams = num2str(lambda_um, '%6.4f');
    pols = ['polaxis'  num2str(optval.polaxis,2)];
    fitswrite(real(fldFull), [mp.full.input_field_rootname '_' lams 'um_' pols '_real.fits' ]);
    fitswrite(imag(fldFull), [mp.full.input_field_rootname '_' lams 'um_' pols '_imag.fits' ]);

    %%--Downsampling for the compact model
    dxF = 1;
    dxC = mp.P1.full.Nbeam/mp.P1.compact.Nbeam;
    Nf = length(fldFull); %--N full
    Nc = ceil_even( (mp.P1.compact.Nbeam/mp.P1.full.Nbeam)*Nf ); %--N compact
    xF = (-Nf/2:Nf/2-1)*dxF;
    xC = (-Nc/2:Nc/2-1)*dxC;
    [Xf, Yf] = meshgrid(xF);
    [Xc, Yc] = meshgrid(xC);
    fldC = interp2(Xf,Yf,fldFull,Xc,Yc,'cubic',0); %--Downsample by interpolation
    fldC = padOrCropEven(fldC,ceil_even(mp.P1.compact.Nbeam+1));
    if(mp.flagPlot)
        figure(607+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
        figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula; drawnow;
    end

    %--Assign to initial E-field in compact model.
    mp.P1.compact.E(:, :, si) = propcustom_relay(fldC, 1, mp.centering);
    
end

% Don't double count the pupil amplitude with the phase retrieval and a model-based mask
mp.P1.compact.mask = ones(size(mp.P1.compact.mask));

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);


%%
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUX RATIO NOISE (FRN) ANALYSIS SECTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Data locations for WFIRST CGI calculations of flux ratio noise (FRN)
mp.path.frn_coro = '~/Downloads/CGdata/'; %--Location of coronagraph performance data tables. Needs slash at end.
fn_prefix = sprintf('s%04dt%04d_',mp.SeriesNum,mp.TrialNum);

%% Change the resolution

E0 = mp.P1.compact.E; %--Don't erase the starting settings.
paths = mp.path;
runLabel = mp.runLabel;
sn = mp.SeriesNum;
tn = mp.TrialNum;
clear mp
mp.runLabel = runLabel;
mp.P1.compact.E = E0; 
mp.path = paths;

%--Re-initialize mp structure
EXAMPLE_defaults_WFIRST_PhaseB_PROPER_SPC_Spec %--Load default model parameters

mp.SeriesNum = sn;
mp.TrialNum = tn;

mp.Fend.res = 5; %--Change the image resolution [pixels per lambda0/D]
mp.full.output_dim = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D

%--Set DM commands back to final
mp.dm1.V = out.dm1.Vall(:,:,end);
mp.dm2.V = out.dm2.Vall(:,:,end);

% %--DEBUGGING:
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.flagParfor = true; %--whether to use parfor for Jacobian calculation

%--Save the config file
fn_config = [mp.path.config mp.runLabel,'_configHD.mat'];
save(fn_config)
fprintf('Saved the config file: \t%s\n',fn_config)
%--Get configuration data from a function file
[mp,out] = falco_init_ws(fn_config);


%% Compute the table of annular zones

mp.eval.Rsens = ...
                [3., 4.;...
                4., 5.;...
                5., 6.;...
                6., 7.;...
                7., 8.]; 
            
tableAnn = falco_FRN_AnnularZone_table(mp);
writetable(tableAnn,[mp.path.frn_coro, fn_prefix, 'AnnZoneList.csv']); %--Save to CSV file
tableAnn  


%% Compute the table InitialRawContrast.csv

[tableContrast, tableCtoNI, data] = falco_FRN_InitialRawContrast(mp);
writetable(tableContrast,[mp.path.frn_coro, fn_prefix, 'InitialRawContrast.csv']); %--Save to CSV file
writetable(tableCtoNI,[mp.path.frn_coro, fn_prefix, 'NItoContrast.csv']); %--Save to CSV file
save([mp.path.frn_coro, fn_prefix, 'c_data.mat'],'data') %--Save 2-D and 1-D Contrast and CtoNI map for making plots later
tableContrast
tableCtoNI

%% Compute the average contrast in each annulus (or annular segment)

Rann = ...
    [3., 4.;...
    4., 5.;...
    5., 8.;...
    8., 9.]; 

Nann = size(Rann,1);

CcohVec = zeros(Nann,1);
CincoVec = zeros(Nann,1);
CtoNIvec = zeros(Nann,1);
rVec = zeros(Nann,1);

for ia=1:Nann        
    min_r = Rann(ia,1);
    max_r = Rann(ia,2);
    rVec(ia) = (min_r+max_r)/2;

    %--Compute the software mask for the scoring region
    maskScore.pixresFP = mp.Fend.res;
    maskScore.rhoInner = min_r; %--lambda0/D
    maskScore.rhoOuter = max_r; %--lambda0/D
    maskScore.angDeg = mp.Fend.score.ang; %--degrees
    maskScore.centering = mp.centering;
    maskScore.FOV = mp.Fend.FOV;
    maskScore.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
    if(isfield(mp.Fend,'shape'));  maskScore.shape = mp.Fend.shape;  end
    [maskPartial,xis,etas] = falco_gen_SW_mask(maskScore);

    %--Compute the average intensity over the selected region
    CcohVec(ia) = sum(sum(maskPartial.*data.Ccoh))/sum(sum(maskPartial));
    CincoVec(ia) = sum(sum(maskPartial.*data.Cinco))/sum(sum(maskPartial));
    CtoNIvec(ia) = sum(sum(maskPartial.*data.CtoNI))/sum(sum(maskPartial));
    % figure(401); imagesc(maskPartial.*Ccoh); axis xy equal tight; colorbar; drawnow; pause(1);
end

%--Outputs for requirements
dataReq.rVec = rVec;
dataReq.CcohVec = CcohVec;
dataReq.CincoVec = CincoVec;
dataReq.CtoNIvec = CtoNIvec;
save([mp.path.frn_coro, fn_prefix, 'c_dataReq.mat'],'dataReq') %--Save 1-D Contrast and CtoNI map for comparing against requirements

%% Compute the Krist table

%--Other constants
mp.yield.Dtel = 2.3631; % meters

%--Define radial sampling and range
mp.yield.R0 = 2.5;
mp.yield.R1 = 9.1;

%--Compute and save the table
tableKrist = falco_FRN_Krist_table(mp);
writetable(tableKrist,[mp.path.frn_coro, fn_prefix, 'KristTable.csv']); %--Save to CSV file

%--Plot the table data
matKrist = tableKrist{:,:};
figure(200); imagesc(matKrist); axis tight;
figure(201); imagesc(log10(matKrist)); axis tight;
figure(202); semilogy(matKrist(:,1),matKrist(:,3),'-b',matKrist(:,1),matKrist(:,4),'-r','Linewidth',3); %--Compare intensity and contrast plots


%% Calculate Sensitivities.csv 

%--Rows 1 to 10: Z2 to Z11 sensitivities to 1nm RMS of Zernike phase aberrations at entrance pupil.
%--Rows 11 to 17: Gain Z5 to Z11 sensitivities
%--Row 18: Pupil X shear
%--Row 19: Pupil Y shear
%--Row 20: DM Settling
%--Row 21: DM Thermal

tableSens = falco_FRN_Sens_table(mp);
writetable(tableSens,[mp.path.frn_coro, fn_prefix, 'Sensitivities.csv']); %--Save to CSV file
tableSens
