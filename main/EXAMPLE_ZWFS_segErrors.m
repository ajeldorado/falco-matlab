% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

clear;

% restoredefaultpath;
% rehash toolboxcache;

%% Step 1: Define Necessary Paths on Your Computer System

pathStem = '/Users/gruane/Desktop/SCDA/';

%--Library locations
mp.path.falco = [pathStem,'falco-matlab/'];  %--Location of FALCO
mp.path.proper = [pathStem,'PROPER/'];  %--Location of FALCO; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = [mp.path.falco,'data/brief/']; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = [mp.path.falco,'data/ws/']; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%---Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

%% Step 2: Load default model parameters

defaults_DST_LC_design

%% Step 3: Overwrite default values as desired

%%--Coronagraph and Pupil Type
mp.coro = 'Roddier';    %--Tested Options: 'LC','HLC','SPLC','Vortex'
mp.flagApod = false;
mp.whichPupil = 'LUVOIR_B_offaxis';

%%--Bandwidth and Wavelength Specs
mp.lambda0 = 400e-9; % central wavelength of bandpass (meters)
mp.fracBW = 10e-9/mp.lambda0 ;%0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 1;  % number of wavelengths or sub-bandpasses (sbp) across entire spectral band 

%%--Pupil definitions

mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.

mp.P1.full.Nbeam = 500; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 

mp.P1.compact.Nbeam = 500;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex.

mp.P4.IDnorm = 0;
mp.P4.ODnorm = 1.1;%0.82;

mp.P1.wGap = 0.01; % Fractional gap width

mp.P4.padFacPct = 0; 


%%-- segmented mirror errors
numSegments = hexSegMirror_numSegments(4); % Number of segments in "full" hex aperture
% LUVOIR B has four rings, but ignores some corner segmentes 

%%-- Focal Plane Mask (F3) Properties

pupilDiam_m = mp.P2.D;%m
maskDepth_m  = 213e-9;%m
maskMaterial = 'FS';%Fused Silica 
maskRadius_lamOverD = 1.06/2;%(maskDiam_m/2)/fNum/mp.lambda0;

mp.F3.Rin = maskRadius_lamOverD;
mp.F3.t = maskDepth_m; % Depth of the Roddier/Zernike spot
mp.FPMampFac = 1.0;% Transmission of the Roddier/Zernike spot
mp.FPMmaterial = maskMaterial;

mp.F3.full.res = 100;
mp.F3.compact.res = 100;


%% Generate the label associated with this trial

mp.runLabel = [mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Get configuration data from a function file

%--Save the config file
fn_config = [mp.path.config mp.runLabel,'_config.mat'];
save(fn_config)
fprintf('Saved the config file: \t%s\n',fn_config)


[mp,out] = falco_init_ws(fn_config);

%% Generate ZWFS images 
disp('***** ZWFS demo *****');
disp('Generating ZWFS calib image...');

b = abs(falco_zwfs_getReferenceWave(mp));
IC0 = falco_zwfs_sim_image(mp);

b = padOrCropEven(b,length(mp.P1.full.mask));
IC0 = padOrCropEven(IC0,length(mp.P1.full.mask));
mask = imerode(logical(mp.P1.full.mask),strel('disk', 2));

%%-- Apply first set of errors
disp('Applying WFE to primary mirror...');
mp.P1.pistons = randn(1,numSegments)/100;% Segment piston in waves 
mp.P1.tiltxs  = randn(1,numSegments)/50;% %Tilts on segments in horiz direction (waves/apDia)
mp.P1.tiltys  = randn(1,numSegments)/50;% %Tilts on segments in vert direction (waves/apDia)

mp = falco_config_gen_chosen_pupil(mp);

actual_phz1 = angle(mp.P1.compact.E(:,:,ceil(mp.Nsbp/2)));

disp('Generating ZWFS image...');
IC1 = falco_zwfs_sim_image(mp);
IC1 = padOrCropEven(IC1,length(mp.P1.full.mask));

%%-- Apply second set of errors
disp('Applying new WFE to primary mirror...');
mp.P1.pistons = mp.P1.pistons + randn(1,numSegments)/2000;% Segment piston in waves 
mp.P1.tiltxs  = mp.P1.tiltxs + randn(1,numSegments)/1000;% %Tilts on segments in horiz direction (waves/apDia)
mp.P1.tiltys  = mp.P1.tiltys + randn(1,numSegments)/1000;% %Tilts on segments in vert direction (waves/apDia)

mp = falco_config_gen_chosen_pupil(mp);

actual_phz2 = angle(mp.P1.compact.E(:,:,ceil(mp.Nsbp/2)));

disp('Generating ZWFS image...');
IC2 = falco_zwfs_sim_image(mp);
IC2 = padOrCropEven(IC2,length(mp.P1.full.mask));

%%-- Reconstruct phases

disp('Reconstructing the wavefront...');
phz1 = falco_zwfs_reconstructor(mp.P1.full.mask,IC1, mask, b,'w');
phz2 = falco_zwfs_reconstructor(mp.P1.full.mask,IC2, mask, b,'w');

phz1 = circshift(rot90(phz1,2),[1 1]);
phz2 = circshift(rot90(phz2,2),[1 1]);
diffphz_meas = phz1 - phz2;

diffphz_actual = actual_phz1 - actual_phz2;
diffphz_actual = padOrCropEven(diffphz_actual,length(diffphz_meas));

%% Plot results

fontsize = 16;

cbmin = -5;
cbmax = 5;

cbmin_diff = -0.2;
cbmax_diff = 0.2;

xvals = -mp.P1.full.Narr/2:mp.P1.full.Narr/2-1;
yvals = xvals';
apRad = mp.P1.full.Nbeam/2;

actual_phz1(actual_phz1==0) = NaN;
actual_phz2(actual_phz2==0) = NaN;

fig0 = figure(666);
set(fig0,'units', 'inches', 'Position', [0 0 12 12])

subplot(3,3,1);
imagesc(xvals/apRad,yvals/apRad,IC0);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Calibration image');

subplot(3,3,2);
imagesc(xvals/apRad,yvals/apRad,IC1);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('WFS image 1');

subplot(3,3,3);
imagesc(xvals/apRad,yvals/apRad,IC2);
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
imagesc(xvals/apRad,yvals/apRad,phz1/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed surf (nm)');

subplot(3,3,8);
imagesc(xvals/apRad,yvals/apRad,phz2/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed surf (nm)');

subplot(3,3,9);
imagesc(xvals/apRad,yvals/apRad,diffphz_actual/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin_diff cbmax_diff]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed difference (nm)');

