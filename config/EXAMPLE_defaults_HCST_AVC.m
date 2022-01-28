
% %--Initialize some structures if they don't already exist

%% Misc

%--Record Keeping
mp.SeriesNum = 867;
mp.TrialNum = 5309;

%--Special Computational Settings
mp.flagParfor = false;
mp.useGPU = false;
mp.flagPlot = false;
mp.flagFiber = false;
mp.flagZWFS = false;

%--General
mp.centering = 'pixel';

%--Whether to include planet in the images
mp.planetFlag = false;

%--Method of computing core throughput:
% - 'HMI' for energy within half-max isophote divided by energy at telescope pupil
% - 'EE' for encircled energy within a radius (mp.thput_radius) divided by energy at telescope pupil
mp.thput_metric = 'EE'; 
mp.thput_radius = 0.7; %--photometric aperture radius [lambda_c/D]. Used ONLY for 'EE' method.
mp.thput_eval_x = 6; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput

%--Where to shift the source to compute the intensity normalization value.
mp.source_x_offset_norm = 7;  % x location [lambda_c/D] in dark hole at which to compute intensity normalization
mp.source_y_offset_norm = 0;  % y location [lambda_c/D] in dark hole at which to compute intensity normalization

%% Bandwidth and Wavelength Specs

mp.lambda0 = 775e-9;    %--Central wavelength of the whole spectral bandpass [meters]
mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing with batch process estimation
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT TESTED YET]
% - 'pwp-iekf' for pairwise probing with iterated extended Kalman filter  [NOT AVAILABLE YET]
mp.estimator = 'perfect';

%--Variables for pairwise probing estimation:
mp.est.flagUseJac = true;    % Whether to use the Jacobian to compute the delta electric fields. If false, the outputs of model_compact are differenced instead.
mp.est.probe.Npairs = 3;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 12;    % Max x/y extent of probed region [actuators].
mp.est.probe.offsetX = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.offsetY = 0;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.

%--New variables for pairwise probing with a Kalman filter
% mp.est.ItrStartKF = 2 %Which correction iteration to start recursive estimate
% mp.est.tExp =
% mp.est.num_im =
% mp.readNoiseStd =
% mp.peakCountsPerPixPerSec =
% mp.est.Qcoef =
% mp.est.Rcoef =
% mp.est.Pcoef0 = 

%% Wavefront Control: General

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppress with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = [];%2:3; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = []; 

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -6:1:1; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
% % mp.ctrl.dm9regfacVec = 1;        %--Additional regularization factor applied to DM9
   
%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1;
mp.dm2.weight = 1;

%--Voltage range restrictions
mp.dm1.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.dm2.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.maxAbsdV = 1000;     %--Max +/- delta voltage step for each actuator for DMs 1 and 2 [volts] %--NOT ENFORCED YET


%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'SM-CVX' for constrained EFC using CVX. --> DEVELOPMENT ONLY
mp.controller = 'gridsearchEFC';

% % % GRID SEARCH EFC DEFAULTS     
%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 20; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1;%1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.dm_ind = [1 2]; %--Which DMs to use


%% Deformable Mirrors: Influence Functions
%--Influence Function Options:
% - 'influence_dm5v2.fits' for one type of Xinetics DM
% - 'influence_BMC_2kDM_400micron_res10.fits' for BMC 2k DM
% - 'influence_BMC_kiloDM_300micron_res10_spline.fits' for BMC kiloDM
mp.dm1.inf_fn = 'influence_BMC_kiloDM_300micron_res10_spline.fits';
mp.dm2.inf_fn = 'influence_BMC_kiloDM_300micron_res10_spline.fits';

mp.dm1.dm_spacing = 300e-6; %--User defined actuator pitch [meters]
mp.dm2.dm_spacing = 300e-6; %--User defined actuator pitch [meters]

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';

%% Deformable Mirrors: Optical Layout Parameters

%--DM1 parameters
mp.dm1.Nact = 34;               % # of actuators across DM array
mp.dm1.VtoH = 1*1e-9*ones(34);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = (34/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = (34/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--DM2 parameters
mp.dm2.Nact = 34;               % # of actuators across DM array
mp.dm2.VtoH = 1*1e-9*ones(34);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;                % clocking of DM surface [degrees]
mp.dm2.xc = (34/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = (34/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm2.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.dm1.Dstop = 100e-3;  %--Diameter of iris [meters]
mp.flagDM2stop = false;  %--Whether to apply an iris or not
mp.dm2.Dstop = 0.4*34e-3;   %--Diameter of iris [meters]

%--DM separations
mp.d_P2_dm1 = 0;        % distance (along +z axis) from P2 pupil to DM1 [meters]
mp.d_dm1_dm2 = 0.20;   % distance between DM1 and DM2 [meters]


%% Optical Layout: All models

%--Key Optical Layout Choices
mp.flagSim = true;      %--Simulation or not
mp.layout = 'Fourier';  %--Which optical layout to use
mp.coro = 'vortex';

%--Final Focal Plane Properties
mp.Fend.res = 5.1; %--Sampling [ pixels per lambda0/D]
mp.Fend.FOV = 12; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin = 3.0;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 120;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = mp.Fend.corr.Rin;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = mp.Fend.corr.Rout;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = mp.Fend.corr.ang;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'left'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 779e-3; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.


%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 256;
mp.P2.compact.Nbeam = 256;
mp.P3.compact.Nbeam = 256;
mp.P4.compact.Nbeam = 256;  % P4 must be the same as P1 for Vortex. 

%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;
mp.NrelayFend = 0; %--How many times to rotate the final image by 180 degrees

% mp.F3.compact.res = 6; % sampling of FPM for compact model [pixels per lambda0/D]

%% Optical Layout: Full Model 

%--Focal Lengths
% mp.fl = 1; 

%--Pupil Plane Resolutions
mp.P1.full.Nbeam = 256;
mp.P2.full.Nbeam = 256;
mp.P3.full.Nbeam = 256;
mp.P4.full.Nbeam = 256;  % P4 must be the same as P1 for Vortex. 

% mp.F3.full.res = 6;    % sampling of FPM for full model [pixels per lambda0/D]

%% Mask Definitions

mp.apodType = 'grayscale';%[mp.path.mask,'ApodizedPupil_500.fits'];
maskaux = fitsread('/Users/jllopsay/Documents/GitHub/falco-matlab/lib/masks/segmentedPupil_noApod.fits');%   ApodizedPupil_HCST
mp.P3.full.mask = imresize(maskaux,[mp.P1.full.Nbeam mp.P1.full.Nbeam]); 
mp.P3.compact.mask = mp.P3.full.mask;
% mp.P3.apodType = 'HCST_AVC';%[mp.path.mask,'ApodizedPupil_500.fits'];

%--Pupil definition
mp.whichPupil = 'Simple';
mp.P1.IDnorm = 0.00; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.ODnorm = 1.00;% Outer diameter of the telescope [diameter]
mp.P1.stretch = 1.00; % factor that stretches the horizontal axis to create elliptical beam 
mp.P1.Dfac = 1; %--Factor scaling inscribed OD to circumscribed OD for the telescope pupil.
mp.P1.Nstrut = 0;% Number of struts 
mp.P1.angStrut = [];%Array of angles of the radial struts (deg)
mp.P1.wStrut = []; % Width of the struts (fraction of pupil diam.)

%--Aperture stop definition
mp.flagApod = true;    %--Whether to use an apodizer or not. Can be a simple aperture stop

%--Lyot stop padding
mp.P4.IDnorm = 0; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.99; %--Lyot stop OD [Dtelescope]
mp.P4.padFacPct = 0;
mp.P4.Nstrut = 0;% Number of struts 
mp.P4.angStrut = [];%Array of angles of the radial struts (deg)
mp.P4.wStrut = []; % Width of the struts (fraction of pupil diam.)

%--Pupil Plane Diameters
mp.P1.D = 15.2e-3/mp.P4.ODnorm; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
mp.P2.D = mp.P1.D;
mp.P3.D = mp.P1.D;
mp.P4.D = 15.2e-3;

%% VC-Specific Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mp.F3.VortexCharge = -8; %--Charge of the vortex mask
%--Pupil definition
mp.whichPupil = 'Simple';
mp.P1.IDnorm = 0; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.ODnorm = 1;% Outer diameter of the telescope [diameter]
mp.P1.Nstrut = 0;% Number of struts 
mp.P1.angStrut = [];%Array of angles of the radial struts (deg)
mp.P1.wStrut = []; % Width of the struts (fraction of pupil diam.)
mp.P1.stretch = 1;

%%-- DM parameters 
mp.dm1.inf_fn = 'influence_BMC_kiloDM_300micron_res10_spline.fits';
mp.dm1.dm_spacing = 3e-4 * 2;  %--User defined actuator pitch
mp.dm1.inf_sign = '+';

xc_best = 16.0335;
yc_best = 16.2271;
xrot_best = 0;
yrot_best = 0;
zrot_best = 179.8671;

mp.dm1.centering = 'pixel';
mp.dm1.Nact = 34;               % # of actuators across DM array
mp.dm1.xtilt = xrot_best;               % for foreshortening. angle of rotation about x-axis [degrees] 
mp.dm1.ytilt = yrot_best;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = zrot_best;           % clocking of DM surface [degrees]
mp.dm1.xc = xc_best;              % x-center location of DM surface [actuator widths]
mp.dm1.yc = yc_best;               % y-center location of DM surface [actuator widths]
Nactbeam = 18.7e-3/mp.dm1.dm_spacing;                % Nactact across the "beam"

% mp.dm1.tied = [441,474];
%mp.dm1.tied = [];

% gains of all actuators [nm/V of free stroke]
% load('/media/hcst/edc9d85f-b1f5-4499-97aa-197291d84a05/HCST_data/efc_results/data/G/gainmap2.mat')
% gainmap = hcst_DM_1Dto2D(bench,gain1tot);
% mp.dm1.VtoH = (4e-7*ones(mp.dm1.Nact) * 2*sqrt(2)).*gainmap';
% mp.aux.VtoH = gainmap;
mp.dm2.inf_fn = 'influence_BMC_kiloDM_300micron_res10_spline.fits';
mp.dm2.dm_spacing = 3e-4 * 2;  %--User defined actuator pitch
mp.dm2.inf_sign = '+';

xc_best = 16.0335;
yc_best = 16.2271;
xrot_best = 0;
yrot_best = 0;
zrot_best = 179.8671;

mp.dm2.centering = 'pixel';
mp.dm2.Nact = 34;               % # of actuators across DM array
mp.dm2.xtilt = xrot_best;               % for foreshortening. angle of rotation about x-axis [degrees] 
mp.dm2.ytilt = yrot_best;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = zrot_best;           % clocking of DM surface [degrees]
mp.dm2.xc = xc_best;              % x-center location of DM surface [actuator widths]
mp.dm2.yc = yc_best;               % y-center location of DM surface [actuator widths]
Nactbeam = 18.7e-3/mp.dm1.dm_spacing;                % Nactact across the "beam"

mp.dm1.HminStep = 0.1;
mp.dm2.HminStep = 0.1;