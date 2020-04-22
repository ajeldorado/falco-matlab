% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

clear;

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command:  addpath(full_path_to_proper); savepath;

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false


%% Step 2: Load default model parameters

EXAMPLE_defaults_DST_LC_design

%% Step 3: Overwrite default values as desired

mp.flagWFS = true; % Activate the WFS mode 
mp.wfs.flagSim = true; % Simulates WFS images, if true 

%%-- Pupil definitions

mp.flagApod = false;
mp.whichPupil = 'LUVOIR_B_offaxis';

mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.

mp.P1.full.Nbeam = 500; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 

mp.P1.compact.Nbeam = 500;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex.

mp.P1.wGap = 0.01; % Fractional gap width
mp.P4.padFacPct = 0; 

%%- segmented mirror errors
numSegments = hexSegMirror_numSegments(4); % Number of segments in "full" hex aperture
% LUVOIR B has four rings, but ignores some corner segments 


%%--ZWFS Mask Properties
mp.wfs.lambda0 = 425e-9;%--Central wavelength of the whole spectral bandpass [meters]
mp.wfs.fracBW = 0.10; %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.wfs.Nsbp = 3; %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control

%%- ZWFS mask properties 
mp.wfs.mask.type = 'transmissive';
maskDepth_m  = 213e-9;%m
maskMaterial = 'FS';%Fused Silica 
maskRadius_lamOverD = 1.06/2;%(maskDiam_m/2)/fNum/mp.lambda0;

mp.wfs.mask.material = maskMaterial; % Required for transmissive mask 
mp.wfs.mask.Rin  = maskRadius_lamOverD;% Radius of the ZWFS dimple 
mp.wfs.mask.Rout = Inf;% Outer 
mp.wfs.mask.depth = maskDepth_m; % Depth of the Zernike dimple
mp.wfs.mask.FPMampFac = 1.0;% Transmission of the Zernike dimple
mp.wfs.mask.res = 20;

%%- ZWFS camera properties 
mp.wfs.cam.Nbeam = mp.P1.full.Nbeam;% beam size at WFS camera 
mp.wfs.cam.Narr = 512; % array size at WFS camera 
mp.wfs.cam.dx = mp.P4.D/mp.wfs.cam.Nbeam;


%% Generate the label associated with this trial

mp.runLabel = [mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Flesh out ws

[mp, out] = falco_flesh_out_workspace(mp);


%% Generate ZWFS images 
rng(1)

disp('***** ZWFS demo *****');
disp('Generating ZWFS calib image...');

Ical = falco_zwfs_getCalibrationImage(mp);
b = falco_zwfs_getReferenceWave(mp);
IZ0 = falco_zwfs_sim_image(mp);

% IZ image mask 
mask = imerode(logical(mp.P1.full.mask),strel('disk', 1));
mask = padOrCropEven(mask,mp.wfs.cam.Narr);


%%-- Apply first set of errors
disp('Applying WFE to primary mirror...');
mp.P1.pistons = randn(1,numSegments)/100;% Segment piston in waves 
mp.P1.tiltxs  = randn(1,numSegments)/50;% %Tilts on segments in horiz direction (waves/apDia)
mp.P1.tiltys  = randn(1,numSegments)/50;% %Tilts on segments in vert direction (waves/apDia)

mp = falco_gen_chosen_pupil(mp);

actual_phz1 = angle(mp.P1.compact.E(:,:,ceil(mp.Nsbp/2)));

disp('Generating ZWFS image...');
IZ1 = falco_zwfs_sim_image(mp);

%%-- Apply second set of errors
disp('Applying new WFE to primary mirror...');
mp.P1.pistons = mp.P1.pistons + randn(1,numSegments)/2000;% Segment piston in waves 
mp.P1.tiltxs  = mp.P1.tiltxs + randn(1,numSegments)/1000;% %Tilts on segments in horiz direction (waves/apDia)
mp.P1.tiltys  = mp.P1.tiltys + randn(1,numSegments)/1000;% %Tilts on segments in vert direction (waves/apDia)

mp = falco_gen_chosen_pupil(mp);

actual_phz2 = angle(mp.P1.compact.E(:,:,ceil(mp.Nsbp/2)));

disp('Generating ZWFS image...');
IZ2 = falco_zwfs_sim_image(mp);

%% Reconstruct phases

disp('Reconstructing the wavefront...');
theta = 2*pi*mp.wfs.mask.depth*(mp.wfs.mask.n(mp.wfs.lambda0)-1)/mp.wfs.lambda0;
phz1 = falco_zwfs_reconstructor(Ical, IZ1, mask, b, theta, 'f');
phz2 = falco_zwfs_reconstructor(Ical, IZ2, mask, b, theta, 'f');

phz1 = circshift(rot90(phz1,2),[1 1]);
phz2 = circshift(rot90(phz2,2),[1 1]);
phz1(~mask) = 0;
phz2(~mask) = 0;

diffphz_meas = phz1 - phz2;
diffphz_meas(~mask) = NaN;

diffphz_actual = actual_phz1 - actual_phz2;
diffphz_actual = padOrCropEven(diffphz_actual,length(diffphz_meas));
diffphz_actual = diffphz_actual - mean(diffphz_actual(mask));
diffphz_actual(~mask) = NaN;

%% Plot results

fontsize = 16;

cbmin = -5;
cbmax = 5;

cbmin_diff = -0.2;
cbmax_diff = 0.2;

xvals = -mp.wfs.cam.Narr/2:mp.wfs.cam.Narr/2-1;
yvals = xvals';
apRad = mp.wfs.cam.Nbeam/2;


fig0 = figure(666);
set(fig0,'units', 'inches', 'Position', [0 0 12 12])

subplot(3,3,1);
imagesc(xvals/apRad,yvals/apRad,IZ0);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Image w/o errors');

subplot(3,3,2);
imagesc(xvals/apRad,yvals/apRad,IZ1);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('WFS image 1');

subplot(3,3,3);
imagesc(xvals/apRad,yvals/apRad,IZ2);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('WFS image 2');

subplot(3,3,4);
imagesc(xvals/apRad,yvals/apRad,actual_phz1/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Actual surf at t1 (nm)');

subplot(3,3,5);
imagesc(xvals/apRad,yvals/apRad,actual_phz2/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Actual surf at t2 (nm)');

subplot(3,3,6);
imagesc(xvals/apRad,yvals/apRad,diffphz_actual/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin_diff cbmax_diff]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Actual difference (nm)');

subplot(3,3,7);
imagesc(xvals/apRad,yvals/apRad,phz1/4/pi*mp.wfs.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed surf (nm)');

subplot(3,3,8);
imagesc(xvals/apRad,yvals/apRad,phz2/4/pi*mp.wfs.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed surf (nm)');

subplot(3,3,9);
imagesc(xvals/apRad,yvals/apRad,diffphz_meas/4/pi*mp.wfs.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin_diff cbmax_diff]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed difference (nm)');

disp(['RMS surface error of difference = ',num2str(rms(diffphz_actual(mask)*mp.lambda0-diffphz_meas(mask)*mp.wfs.lambda0)/4/pi*1e12),' pm']);