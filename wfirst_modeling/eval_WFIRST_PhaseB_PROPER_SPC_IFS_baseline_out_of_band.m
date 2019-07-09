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

EXAMPLE_defaults_WFIRST_PhaseB_PROPER_SPC_IFS


%% Step 3: Overwrite default values as desired

% %%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 45;
mp.TrialNum = 1;


%% Change the wavelength and resolution

%--Testing
mp.lambda0 = 730e-9;   %--Central wavelength of the whole spectral bandpass [meters]
mp.fracBW = 0.15;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 7;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

% %--Band 1
% mp.lambda0 = 575e-9;   %--Central wavelength of the whole spectral bandpass [meters]
% mp.fracBW = 0.033;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 5;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
% 
% %--H-alpha
% mp.lambda0 = 656e-9;   %--Central wavelength of the whole spectral bandpass [meters]
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 3;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
% 
% %--Band 4
% mp.lambda0 = 825e-9;   %--Central wavelength of the whole spectral bandpass [meters]
% mp.fracBW = 0.033;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 5;          %--Number of wavelengths to used to approximate an image in each sub-bandpass



mp.full.lambda0_m = mp.lambda0;
mas2lam0D = 1/(mp.lambda0/2.3631*180/pi*3600*1000);
mp.Fend.res = 1/(21*mas2lam0D); %--Change the image resolution [pixels per lambda0/D]. Detector resolution is 21 mas/pixel.
mp.Fend.FOV = 15.; %--half-width of the field of view in both dimensions [lambda0/D]

mp.full.output_dim = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D

%--Set DM commands back to the dark hole settings
load('Series0045_Trial0001_SPLC_WFIRST180718_2DM48_z1_IWA2.6_OWA9_5lams730nm_BW15_plannedEFC_snippet.mat','out');
mp.dm1.V = out.dm1.Vall(:,:,end);
mp.dm2.V = out.dm2.Vall(:,:,end);

%% Add a sinusoid to DM1
% 
% xs = (1:48)/(46.3);
% [XS,YS] = meshgrid(xs);
% 
% w = 6; %--offset
% 
% sinusoid = sin(2*pi*XS*w*cosd(mp.dm1.ytilt));
% 
% figure(250); imagesc(sinusoid); axis xy equal tight; colorbar;
% 
% c = 1e-7;
% dV = sqrt(c)*(4*pi*mp.lambda0)./mp.dm1.VtoH.*sinusoid;
% mp.dm1.V = mp.dm1.V + dV;
% 
% figure(251); imagesc(dV); axis xy equal tight; colorbar; set(gca,'Fontsize',20);


%% Step 3b: Obtain the phase retrieval phase.

mp.full.input_field_rootname = '/home/ajriggs/Repos/falco-matlab/data/maps/input_full';
optval = mp.full;
optval.source_x_offset =0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.dm1_m = fitsread([mp.full.data_dir 'errors_polaxis10_dm.fits']);
optval.use_dm1 = 1;

optval.end_at_fpm_exit_pupil = 1;
optval.output_field_rootname = [fileparts(mp.full.input_field_rootname) filesep 'fld_at_xtPup'];
optval.use_fpm = 0;
optval.use_hlc_dm_patterns = 0;
nout = 1024; %512; 			% nout > pupil_daim_pix
optval.output_dim = 1024;%% Get the Input Pupil's E-field

if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

%--Get the Input Pupil's E-field
mp.P1.compact.E = ones(mp.P1.compact.Nbeam+2,mp.P1.compact.Nbeam+2,mp.Nsbp); %--Initialize
for si=1:mp.Nsbp
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
    [Xf,Yf] = meshgrid(xF);
    [Xc,Yc] = meshgrid(xC);
    fldC = interp2(Xf,Yf,fldFull,Xc,Yc,'cubic',0); %--Downsample by interpolation
    fldC = padOrCropEven(fldC,ceil_even(mp.P1.compact.Nbeam+1));
    if(mp.flagPlot)
        figure(607+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
        figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula; drawnow;
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

%% Set up the rest of the workspace

%--Save the config file
fn_config = [mp.path.config mp.runLabel,'_config_' num2str(round(mp.lambda0*1e9)) 'nm_BW' num2str(round(100*mp.fracBW)) '.mat'];
save(fn_config)
fprintf('Saved the config file: \t%s\n',fn_config)
%--Get configuration data from a function file
[mp,~] = falco_init_ws(fn_config);


%% Take broadband image 

Im = falco_get_summed_image(mp);

%%

figure(900)
imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im),[-9 -6]); 
axis xy equal tight; ch_psf=colorbar; colormap parula;
xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
ylabel(ch_psf,'$log_{10}$(NI)','Fontsize',24,'Interpreter','LaTex');
title(sprintf('PSF for %d%% BW at %d nm',round(100*mp.fracBW),round(1e9*mp.lambda0)),'Fontsize',20,'Fontweight','Bold');
set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
set(gcf,'Color','w')

fn = sprintf('/home/ajriggs/Downloads/s45t01/PSF_SPC-IFS_BW%02dat%dnm.fits',round(100*mp.fracBW),round(1e9*mp.lambda0));
fitswrite(Im,fn);


