% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to:
%  1) Specify the key parameter values for a Lyot coronagraph.
%  2) Load the rest of the default settings.
%  3) Save out all the input parameters.
%  4) Run a single trial of WFSC using FALCO.


clear;

restoredefaultpath;
rehash toolboxcache;

%% Define Necessary Paths on Your System

pathStem = '/Users/gruane/Desktop/JPL/SCDA/';

%--Library locations
mp.path.falco = [pathStem,'falco-matlab/'];  %--Location of FALCO
mp.path.proper = [pathStem,'PROPER/'];  %--Location of FALCO; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = [mp.path.falco,'data/brief/']; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = [mp.path.falco,'data/ws/']; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


%% Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

mp.flagPlot = true;

%%

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing with batch process estimation
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT AVAILABLE YET]
% - 'pwp-iekf' for pairwise probing with iterated extended Kalman filter  [NOT AVAILABLE YET]
% mp.estimator = 'pwp-bp';
mp.controller = 'perfect';

%--New variables for estimation:
% - Note: For 360-degree dark hole, must set mp.est.probe.Npairs>=3 and mp.est.probe.axis = 'alternate'.
mp.est.probe.Npairs = 3;%2;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 50;%12;%20;    % Max x/y extent of probed region [actuators].
mp.est.probe.offsetX = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.offsetY = 0;   % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';    % which axis to have the phase discontinuity along [x or y or xy/alt/alternate].
mp.est.probe.gainFudge = 6;%4;%1;     % empirical fudge factor to make average probe amplitude match desired value.

mp.flagSim = true;
mp.layout = 'Fourier';

% mp.dm1.inf_fn = 'influence_BMC_kiloDM_300um_N131.fits';
% mp.dm2.inf_fn = 'influence_BMC_kiloDM_300um_N131.fits';
mp.dm1.inf_fn = 'influence_BMC_kiloDM_300um_N65.fits';
mp.dm2.inf_fn = 'influence_BMC_kiloDM_300um_N65.fits';
% mp.dm1.inf_fn = 'influence_dm5v2.fits';
% mp.dm2.inf_fn = 'influence_dm5v2.fits';

mp.dm1.dm_spacing = 1e-3; %--User defined actuator pitch
mp.dm2.dm_spacing = 1e-3; %--User defined actuator pitch

mp.P2.D =     30e-3; % beam diameter at pupil closest to the DMs  (meters)
mp.dm1.Nact = 32; % number of actuators across DM1
mp.dm2.Nact = 32; % number of actuators across DM2
mp.dm_weights = ones(9,1);   % vector of relative weighting of DMs for EFC

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';

%% Step 1: Define any variable values that will overwrite the defaults (in falco_config_defaults_SPLC)

% %%--Record Keeping
% mp.TrialNum = 1; %--Always use a diffrent Trial # for different calls of FALCO.
% mp.SeriesNum = 1; %--Use the same Series # for sets of similar trials.

mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

% % Controller options: 
% %  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
% %  - 'plannedEFC' for EFC with an automated regularization schedule
% %  - 'conEFC' for constrained EFC using CVX. --> DEVELOPMENT ONLY
% mp.controller = 'gridsearchEFC'; %'plannedEFC';%--Controller options: 'gridsearchEFC' or 'plannedEFC'

% % Estimator options:
% mp.estimator = 'perfect'; 
% 

%%--Coronagraph and Pupil Type
mp.coro = 'Roddier';    %--Tested Options: 'LC','HLC','SPLC','Vortex'
mp.flagApod = false;
mp.whichPupil = 'LUVOIR_B_offaxis';

%--Zernikes to suppress with controller
mp.jac.zerns = 1;%[1 2 3]; %1; %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include at least 1 for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = 2:3; %--Noll indices of Zernikes to compute values for
mp.eval.Rsens = [3, 4;... %--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
                 4, 8];    


%%--Pupil Plane and DM Plane Properties
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)
mp.d_dm1_dm2 = 3; % distance between DM1 and DM2 (meters)


%%--Cases for different bandwidth specifications:
% Case 1: Design,Modeling,Testbed. 1 wvl per bandpass. Same wvls in compact and full models
  %--> For testbed, make wvl centers the centers of the sbp instead of
  %linearly from end to end. 
% Case 2: Modeling/Testbed w/ estimation. 1 wvl/sbp in compact model. Nwpsbp wvl/sbp in full model. Need to treat end wavelengths specially in full model.

%%--Bandwidth and Wavelength Specs
% NOTE: Actual wavelengths can be set directly in "mp.sbp_centers" vector in 
% falco_init_ws.m instead if desired. Otherwise, they are auto-generated 
% from the bandwidth and # of wavelengths specified below.
mp.lambda0 = 400e-9; % central wavelength of bandpass (meters)
mp.fracBW = 10e-9/mp.lambda0 ;%0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 1;  % number of wavelengths or sub-bandpasses (sbp) across entire spectral band 

% To map the Lyot stop and DM
f_OX6 = 357.394e-3;
f_OCshrt = 399.338e-3;
f_OClong = 799.032e-3;
magn = f_OX6/f_OClong;


%% -Pupil definitions

mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.

mp.P1.full.Nbeam = 500; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 

mp.P1.compact.Nbeam = 500;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex.

% mp.P3.apodType = 'Simple';
% 
% mp.P3.IDnorm = 0;
% mp.P3.ODnorm = 0.84;

mp.P4.IDnorm = 0;
mp.P4.ODnorm = 1.1;%0.82;

%mp.P1.gapWidth = 6e-3/mp.P1.D; % Fractional gap width
mp.P1.gapWidth = 0.01; % Fractional gap width

mp.P4.padFacPct = 0; 


%% segmented mirror errors
numSegments = hexSegMirror_numSegments(4); % Number of segments in "full" hex aperture
% LUVOIR B has four rings, but ignores some corner segmentes 

% Note: the following are defined for lambda0
% mp.P1.pistons = randn(1,numSegments)/100;% Segment piston in waves 
% mp.P1.tiltxs  = randn(1,numSegments)/50;% %Tilts on segments in horiz direction (waves/apDia)
% mp.P1.tiltys  = randn(1,numSegments)/50;% %Tilts on segments in vert direction (waves/apDia)

%% Focal Plane Mask (F3) Properties

pupilDiam_m = mp.P2.D;%m
% maskDiam_m = 19.5e-6;%m
maskDepth_m  = 213e-9;%m
maskMaterial = 'FS';%Fused Silica 
% fNum = f_OCshrt/pupilDiam_m;
maskRadius_lamOverD = 1.06/2;%(maskDiam_m/2)/fNum/mp.lambda0;

mp.F3.Rin = maskRadius_lamOverD;
mp.F3.t = maskDepth_m; % Depth of the Roddier/Zernike spot
mp.FPMampFac = 1.0;% Transmission of the Roddier/Zernike spot
mp.FPMmaterial = maskMaterial;

mp.F3.full.res = 100;
mp.F3.compact.res = 100;


%% Final Focal Plane (Fend. Properties


%--Specs for Correction (Corr) region and the Scoring (Score) region.
mp.Fend.corr.Rin  = 4; %--lambda0/D, inner radius of correction region
mp.Fend.score.Rin = mp.Fend.corr.Rin; %--Needs to be <= that of Correction mask
mp.Fend.corr.Rout  = 10;%floor(mp.dm1.Nact/2*(1-mp.fracBW/2)); %--lambda0/D, outer radius of correction region
mp.Fend.score.Rout = mp.Fend.corr.Rout; %--Needs to be <= that of Correction mask
mp.Fend.corr.ang  = 180; %--degrees per side
mp.Fend.score.ang = 180; %--degrees per side
mp.Fend.sides = 'right'; %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


% %%--Final Focal Plane (Fend. Properties
mp.Fend.res = 4; %--Pixels per lambda_c/D
% mp.Fend.FOV = 1 + mp.Fend.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D

%% Part 2: Call the function to define the rest of the variables and initialize the workspace
if(exist('mp','var')==false); mp.dummy = 1; end %--Initialize the structure mp if it doesn't exist already.

mp = falco_config_defaults_Roddier(mp); %--Load defaults for undefined values

%% Part 3: Generate the label associated with this trial

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

b = abs(falco_zwfs_getReferenceWave(mp));
IC0 = falco_zwfs_sim_image(mp);

b = padOrCropEven(b,length(mp.P1.full.mask));
IC0 = padOrCropEven(IC0,length(mp.P1.full.mask));
mask = imerode(logical(mp.P1.full.mask),strel('disk', 2));

%%-- Apply first set of errors
mp.P1.pistons = randn(1,numSegments)/100;% Segment piston in waves 
mp.P1.tiltxs  = randn(1,numSegments)/50;% %Tilts on segments in horiz direction (waves/apDia)
mp.P1.tiltys  = randn(1,numSegments)/50;% %Tilts on segments in vert direction (waves/apDia)

mp = falco_config_gen_chosen_pupil(mp);

actual_phz1 = angle(mp.P1.compact.E(:,:,ceil(mp.Nsbp/2)));

IC1 = falco_zwfs_sim_image(mp);
IC1 = padOrCropEven(IC1,length(mp.P1.full.mask));

%%-- Apply second set of errors
mp.P1.pistons = mp.P1.pistons + randn(1,numSegments)/2000;% Segment piston in waves 
mp.P1.tiltxs  = mp.P1.tiltxs + randn(1,numSegments)/1000;% %Tilts on segments in horiz direction (waves/apDia)
mp.P1.tiltys  = mp.P1.tiltys + randn(1,numSegments)/1000;% %Tilts on segments in vert direction (waves/apDia)

mp = falco_config_gen_chosen_pupil(mp);

actual_phz2 = angle(mp.P1.compact.E(:,:,ceil(mp.Nsbp/2)));

IC2 = falco_zwfs_sim_image(mp);
IC2 = padOrCropEven(IC2,length(mp.P1.full.mask));

%%-- Reconstruct phases
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
title('Actual surf at t1');

subplot(3,3,5);
imagesc(xvals/apRad,yvals/apRad,actual_phz2/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Actual surf at t2');

subplot(3,3,6);
imagesc(xvals/apRad,yvals/apRad,diffphz_actual/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin_diff cbmax_diff]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Actual difference');

subplot(3,3,7);
imagesc(xvals/apRad,yvals/apRad,phz1/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed surf at t1');

subplot(3,3,8);
imagesc(xvals/apRad,yvals/apRad,phz2/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin cbmax]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed surf at t2');

subplot(3,3,9);
imagesc(xvals/apRad,yvals/apRad,diffphz_actual/4/pi*mp.lambda0*1e9);
colormap(parula(256));hcb=colorbar;
axis image;%axis off; 
axis([-1 1 -1 1]);caxis([cbmin_diff cbmax_diff]);
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1,'FontSize', fontsize);
set(gca,'TickDir','out');set(gca,'YDir','normal');
title('Reconstructed difference');

