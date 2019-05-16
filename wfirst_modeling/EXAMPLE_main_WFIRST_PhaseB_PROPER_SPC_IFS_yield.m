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

mp.Fend.res = 5;
mp.full.output_dim = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D

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
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation

% %--DEBUGGING:
% mp.fracBW = 0.04;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 2;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.flagParfor = true; %--whether to use parfor for Jacobian calculation

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

mp.P1.compact.E = zeros(mp.P1.compact.Nbeam+2,mp.P1.compact.Nbeam+2,mp.Nsbp); %--Initialize


%--DEBUGGING: Perfect wavefront
mp.full.pol_conds = [0];
mp.full.polaxis = 0;
mp.full.use_errors = false;
mp.full.dm1.flatmap = 0;

optval.phaseb_dir = mp.full.phaseb_dir;

optval.cor_type = mp.full.cor_type;

optval.source_x_offset =0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.use_errors = mp.full.use_errors;
optval.polaxis = mp.full.polaxis; %-2; 
% 
% % 2. full model, for regular psf 
% % EE = prop_run_multi(['wfirst_phaseb_v2'], lam_array, npsf, 'quiet', 'passvalue',optval );
% EE = prop_run(['wfirst_phaseb_v2'], lam_array, npsf, 'quiet', 'passvalue',optval );

% 3. full model, for field

optval.dm1_m = mp.full.dm1.flatmap; %fitsread([mp.full.phaseb_dir 'dm1_flatten.fits']);
optval.use_dm1 =1;

optval.end_at_fpm_exit_pupil = 1;
optval.output_field_rootname = ['fld_at_xtPup'];
optval.use_fpm=0;
optval.use_hlc_dm_patterns=0;
nout = 1024;%512; 			% nout > pupil_daim_pix
% if testcase >2; nout =1024; end % >= pupil_daim_pix; 
% fld = prop_run_multi(['wfirst_phaseb_v2'], lam_array, nout, 'quiet', 'passvalue',optval );

%% Check different polarization state initializations.
% mp.full.pol_conds = [-2,-1,1,2];
% 
% for ipol=1:length(mp.full.pol_conds)
%     lambda_um = 1e6*mp.lambda0;%*lambdaFacs(si);
%     
%     optval.polaxis = mp.full.pol_conds(ipol);
% 
%     fld = prop_run(['wfirst_phaseb_v2b'], lambda_um, nout, 'quiet', 'passvalue',optval );
%     % % % fld(2:end,2:end) = rot90(fld(2:end,2:end),2);
% 
%     % figure(601); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv;
%     % figure(602); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula;
%     figure(500+ipol); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv; drawnow;
%     figure(600+ipol); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula; drawnow;
%     
% end
% %%
% return
%%

if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

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

    %--Subtract out the tip/tilt
    Narray = ceil_even(mp.P1.compact.Nbeam+1);

    x = (-Narray/2:Narray/2-1)/mp.P1.compact.Nbeam; 
    [tip,tilt] = meshgrid(x);
    
    fldC = padOrCropEven(fldC,ceil_even(mp.P1.compact.Nbeam+1));
    if(mp.full.dm1.flatmap ~= 0)
        fldC = fldC.*exp(1i*tilt*0.87*(1/lambdaFacs(si)));
    end
    
    temp = 0*fldC;
    temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
    mp.P1.compact.E(:,:,si) = temp;
    
    figure(617+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
    
    

end

% %% Find the coefficient of the tip to use on the phase retrieval map to flatten it
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

[out,mp] = falco_wfsc_loop(mp);

%% Initialize: Sort out file paths and save the config file    

% %--Add the slash or backslash to the FALCO path if it isn't there.
% if( strcmp(mp.path.falco(end),'/')==false && strcmp(mp.path.falco(end),'\')==false )
%     mp.path.falco = [mp.path.falco filesep];
% end
% 
% mp.path.dummy = 1; %--Initialize the folders structure in case it doesn't already exist
% 
% %--Store minimal data to re-construct the data from the run: the config files and "out" structure after a trial go here
% if(isfield(mp.path,'config')==false)
%     mp.path.config = [mp.path.falco filesep 'data' filesep 'config' filesep];     
% end
% 
% %--Entire final workspace from FALCO gets saved here.
% if(isfield(mp.path,'ws')==false)
%     mp.path.ws = [mp.path.falco filesep 'data' filesep 'ws' filesep];      
% end
% 
% %--Save the config file
% fn_config = [mp.path.config mp.runLabel,'_config.mat'];
% save(fn_config)
% fprintf('Saved the config file: \t%s\n',fn_config)
% 
% %%--Get configuration data from a function file
% % if(~mp.flagSim);  bench = mp.bench;  end %--Save the testbed structure "mp.bench" into "bench" so it isn't overwritten by falco_init_ws
% [mp,out] = falco_init_ws(fn_config);
% % if(~mp.flagSim);  mp.bench = bench;  end






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Change the resolution

E0 = mp.P1.compact.E; %--Don't erase the starting settings.
paths = mp.path;
runLabel = mp.runLabel;
clear mp
mp.runLabel = runLabel;
mp.P1.compact.E = E0; 
mp.path = paths;

%--Re-initialize mp structure
EXAMPLE_defaults_WFIRST_PhaseB_PROPER_SPC_IFS %--Load default model parameters

mp.Fend.res = 5; %--Change the image resolution [pixels per lambda0/D]
mp.full.output_dim = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D

%--Set DM commands back to final
mp.dm1.V = out.dm1.Vall(:,:,end);
mp.dm2.V = out.dm2.Vall(:,:,end);

%--DEBUGGING:
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation


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
            
%--Annular Zone Table
Nann = size(mp.eval.Rsens,1);
tableAnn = zeros(2*Nann,3);
tableAnn(:,1) = [(0:(Nann-1)).';(0:(Nann-1)).']; %--First column
tableAnn(Nann+1:end,2) = 1; %--Second column
tableAnn(1:Nann,3) = mp.eval.Rsens(:,1); %--Third column, 1st half
tableAnn(Nann+1:end,3) = mp.eval.Rsens(:,2); %--Third column, 2nd half
tableAnn

% [tableAnn,fnAnn] = falco_gen_FRNtable_AnnZoneList(mp);

%% Compute the table InitialRawContrast.csv --> DO THIS INSIDE OF THE FRN CALCULATOR TO RE-USE THE CONTRAST MAPS

            
% [tableRaw,fnRaw] = falco_gen_FRNtable_InitialRawContrast(mp);

%--Initial Raw Contrast
Nann = size(mp.eval.Rsens,1);
tableRaw = zeros(4*Nann,5);
tableRaw(:,1) = 0:(4*Nann-1); %--First column
tableRaw(2:2:end,2) = 1; %--Second column
tableRaw(1:2:end,3) = [(0:(Nann-1)).';(0:(Nann-1)).']; %--Third column, one half
tableRaw(2:2:end,3) = [(0:(Nann-1)).';(0:(Nann-1)).']; %--Third column, other half
tableRaw(2*Nann+1:end,4) = 1; %--Fourth Column
tableRaw

%--Compute coherent-only image
mp.full.pol_conds = 10; %--Which polarization states to use when creating an image.
ImCoh = falco_get_summed_image(mp);

%--Compute coherent+incoherent light image
mp.full.pol_conds = [-2,-1,1,2]; %--Which polarization states to use when creating an image.
ImBoth = falco_get_summed_image(mp);
ImInco = ImBoth-ImCoh;

%% Compute 2-D intensity-to-contrast map. 
%  Compute at the center wavelength with the compact model

%--Compute the (x,y) pairs in the image for each pixel
[XIS,ETAS] = meshgrid(mp.Fend.xisDL,mp.Fend.etasDL);

maskBoolQuad1 = mp.Fend.corr.maskBool & XIS>=0 & ETAS>=0;
figure(89); imagesc(maskBoolQuad1); axis xy equal tight; colorbar;

coords = [XIS(maskBoolQuad1),ETAS(maskBoolQuad1)];
Nc = size(coords,1);

tic
peakVals = zeros(Nc,1);
if(mp.flagParfor)
    parfor ic = 1:Nc
        peakVals(ic) = falco_get_offset_peak(mp,coords,ic);
    end
else
   for ic = 1:Nc
        modvar.whichSource = 'offaxis';
        modvar.x_offset = coords(ic,1); 
        modvar.y_offset = coords(ic,2); 
        modvar.sbpIndex = mp.si_ref; 
        modvar.wpsbpIndex = mp.wi_ref;
        E2D = model_compact(mp, modvar);
        peakVals(ic) = max(abs(E2D(:)).^2);
    end  
end
fprintf('Time = %.2f s\n',toc)

peak2D = zeros(size(ETAS));
peak2D(maskBoolQuad1) = peakVals;

peak2D(2:mp.Fend.Neta/2,:) = flipud(peak2D(mp.Fend.Neta/2+2:end,:)); %--Fill in Quadrant 4
peak2D(:,2:mp.Fend.Neta/2) = fliplr(peak2D(:,mp.Fend.Neta/2+2:end)); %--Fill in Quadrants 2 and 3

CtoNI = peak2D/max(peakVals);

figure(90); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,peak2D); axis xy equal tight; colorbar;
figure(91); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,CtoNI); axis xy equal tight; colorbar;

%%
%--Compute the average contrast in each sector or annulus
% int_vec = zeros(Noff,1);
% for ioff=1:Noff        
%     min_r = rs(ioff) - dr/2;
%     max_r = rs(ioff) + dr/2;
% 
%     %--Compute the software mask for the scoring region
%     maskScore.pixresFP = mp.Fend.res;
%     maskScore.rhoInner = min_r; %--lambda0/D
%     maskScore.rhoOuter = max_r; %--lambda0/D
%     maskScore.angDeg = mp.Fend.score.ang; %--degrees
%     maskScore.centering = mp.centering;
%     maskScore.FOV = mp.Fend.FOV;
%     maskScore.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
%     if(isfield(mp.Fend,'shape'));  maskScore.shape = mp.Fend.shape;  end
%     [maskPartial,xis,etas] = falco_gen_SW_mask(maskScore);
% 
%     %--Compute the average intensity over the selected region
%     int_vec(ioff) = sum(sum(maskPartial.*Icam))/sum(sum(maskPartial));
% end

%% Compute the Krist table


% %--Other constants
% % lambda_center = mp.lambda0;
% mp.yield.Dtel = 2.3631; % meters
% 
% %--Define radial sampling and range
% % pixel_step = 1/mp.Fend.eval.res; %1/mp.Fend.res; % step size between pixels[lambda0/D]
% mp.yield.R0 = 2.5;
% mp.yield.R1 = 9.1;
% 
%  
% tableOut = falco_compute_FRN_input_table(mp);
% 
% figure(200); imagesc(tableOut); axis tight;
% figure(201); imagesc(log10(tableOut)); axis tight;

%% Calculate sensitivities to 1nm RMS of Zernike phase aberrations at entrance pupil.

mp.full.ZrmsVal = 1e-9; %--RMS values for each Zernike specified in vector indsZnoll [meters] 
mp.eval.Rsens = ...
                [3., 4.;...
                4., 5.;...
                5., 6.;...
                6., 7.;...
                7., 8.]; 
mp.eval.indsZnoll = 2:3; %[2,3];

if( isempty(mp.eval.Rsens)==false || isempty(mp.eval.indsZnoll)==false )
    Zsens = falco_get_Zernike_sensitivities(mp);
end

%%


function peakVal = falco_get_offset_peak(mp,coords,ic)
    
    modvar.whichSource = 'offaxis';
    modvar.x_offset = coords(ic,1); 
    modvar.y_offset = coords(ic,2); 
    modvar.sbpIndex = mp.si_ref; 
    modvar.wpsbpIndex = mp.wi_ref;
    E2D = model_compact(mp, modvar);
    peakVal = max(abs(E2D(:)).^2);
    
end

