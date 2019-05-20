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



%% Step 3a) Define custom SPC values.

%--Bandwidth and Wavelength Specs
mp.lambda0 = 730e-9;   %--Central wavelength of the whole spectral bandpass [meters]
mp.Nwpsbp = 3;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%--Full model mask files and resolutions
mp.full.cor_type = 'spc_ifs_custom';
file_dir = '/home/ajriggs/Documents/Sim/cgi/wfirst_phaseb/spc_ifs_custom/';
mp.full.pupil_mask_file = [file_dir, 'SPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW15Nlam6.fits'];        mp.fracBW = 0.15; mp.Nsbp = 5;%--SPM file name
% mp.full.pupil_mask_file = [file_dir, 'minpadSPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW15Nlam6.fits'];  mp.fracBW = 0.15; mp.Nsbp = 5;%--SPM file name
% mp.full.pupil_mask_file = [file_dir, 'SPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW2Nlam6.fits'];         mp.fracBW = 0.02; mp.Nsbp = 1;%--SPM file name
% mp.full.pupil_mask_file = [file_dir, 'minpadSPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW2Nlam6.fits'];   mp.fracBW = 0.02; mp.Nsbp = 1;%--SPM file name

mp.full.pupil_file = [file_dir, 'pupil_SPC-20190130_rotated.fits'];
mp.full.pupil_diam_pix = 1000;
mp.full.fpm_file = [file_dir, 'fpm_sharp_26WA90_65deg_res20.fits'];
mp.full.fpm_sampling_lam0 = 0.05; 	% sampling in lambda0/D of FPM mask
mp.full.lyot_stop_file = [file_dir,'LS_sharp_20D96_95deg_N500.fits'];
mp.full.lambda0_m = mp.lambda0;       % FPM scaled for this central wavelength

%--Compact model FPM
mp.compact.flagGenFPM = true;
mp.F3.Rin = 2.6;   % inner hard-edge radius of the focal plane mask [lambda0/D]. Needs to be <= mp.F3.Rin 
mp.F3.Rout = 9;   % radius of outer opaque edge of FPM [lambda0/D]
mp.F3.ang = 65;    % on each side, opening angle [degrees]
mp.F3.compact.mask = rmfield(mp.F3.compact.mask,'amp');
mp.F3.compact = rmfield(mp.F3.compact,'mask');

%--Lyot stop shape
mp.compact.flagGenLS = true;
mp.LSshape = 'bowtie';
mp.P4.IDnorm = 0.20; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.96; %--Lyot stop OD [Dtelescope]
mp.P4.ang = 95;      %--Lyot stop opening angle [degrees]
mp.P4.wStrut = 0;    %--Lyot stop strut width [pupil diameters]


%% Step 3b: Overwrite default values as desired

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
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% % mp.flagParfor = false; %--whether to use parfor for Jacobian calculation


mp.ctrl.sched_mat = [...
    [0,0,0,1,0];
    repmat([1,1j,12,0,1],[5,1]);...
    [1,-5,12,0,0];...
    repmat([1,1j,12,0,1],[10,1]);...
    ];
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

% %--GRID SEARCH EFC    
% mp.controller = 'gridsearchEFC';
% mp.Nitr = 5; %--Number of estimation+control iterations to perform
% mp.relinItrVec = 1;  %--Which correction iterations at which to re-compute the control Jacobian
% mp.ctrl.flagUseModel = true; %--Whether to perform a model-based (vs empirical) grid search for the controller

%%
%--Shaped Pupil Mask: Load and downsample.
mp.SPname = 'SPC-custom';
SP0 = fitsread(mp.full.pupil_mask_file);
% % SP0(2:end,2:end) = rot90(SP0(2:end,2:end),2);

if(mp.P1.compact.Nbeam==1000)
    mp.P3.compact.mask = SP0;
else
    SP0 = SP0(2:end,2:end);
    % figure(1); imagesc(SP0); axis xy equal tight; colormap jet; colorbar;
    % figure(11); imagesc(SP0-fliplr(SP0)); axis xy equal tight; colormap jet; colorbar;
    dx0 = 1/1000;
    dx1 = 1/mp.P3.compact.Nbeam;
    N0 = size(SP0,1);
    switch lower(mp.centering)
        case{'pixel'}
            N1 = ceil_odd(N0*dx0/dx1);
        case{'interpixel'}
            N1 = ceil_even(N0*dx0/dx1);
    end
    x0 = (-(N0-1)/2:(N0-1)/2)*dx0;
    [X0,Y0] = meshgrid(x0);
    R0 = sqrt(X0.^2+Y0.^2);
    Window = 0*R0;
    Window(R0<=dx1) = 1; Window = Window/sum(sum(Window));
    % figure(10); imagesc(Window); axis xy equal tight; colormap jet; colorbar;
    SP0 = ifftshift(  ifft2( fft2(fftshift(Window)).*fft2(fftshift(SP0)) )); %--To get good grayscale edges, convolve with the correct window before downsampling.
    SP0 = circshift(SP0,[1 1]); %--Undo a centering shift
    x1 = (-(N1-1)/2:(N1-1)/2)*dx1;
    [X1,Y1] = meshgrid(x1);
    SP1 = interp2(X0,Y0,SP0,X1,Y1,'cubic',0); %--Downsample by interpolation

    switch lower(mp.centering)
        case{'pixel'}
            mp.P3.compact.mask = zeros(N1+1,N1+1);
            mp.P3.compact.mask(2:end,2:end) = SP1;
        otherwise
            mp.P3.compact.mask = SP1;
    end
    
    figure(2); imagesc(SP0); axis xy equal tight; colormap jet; colorbar; drawnow;
    figure(3); imagesc(SP1); axis xy equal tight; colormap jet; colorbar; drawnow;
    % figure(12); imagesc(SP0-fliplr(SP0)); axis xy equal tight; colormap jet; colorbar;
    % figure(13); imagesc(SP1-fliplr(SP1)); axis xy equal tight; colormap jet; colorbar;

end


%% Step 3c: Obtain the phase retrieval phase.

mp.full.input_field_rootname = '/home/ajriggs/Repos/falco-matlab/data/maps/input_full';


optval = mp.full;

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
optval.use_fpm = 0;
optval.use_hlc_dm_patterns = 0;
nout = 1024;%512; 			% nout > pupil_daim_pix

optval.output_dim = 1024;

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

%     x = (-Narray/2:Narray/2-1)/mp.P1.compact.Nbeam; 
%     [tip,tilt] = meshgrid(x);
    
    fldC = padOrCropEven(fldC,ceil_even(mp.P1.compact.Nbeam+1));
%     fldC = fldC.*exp(1i*tilt*0.87*(1/lambdaFacs(si)));
    
    temp = 0*fldC;
    temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
    mp.P1.compact.E(:,:,si) = temp;
    
    figure(617+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
    
end

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);
