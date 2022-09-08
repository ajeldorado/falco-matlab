% Initialize some structures if they don't already exist.
% Generate or load the pupil masks.

%% Misc

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;
mp.runLabel = sprintf('Series%04d_Trial%04d', mp.SeriesNum, mp.TrialNum);

%--Special Computational Settings
mp.flagParfor = false;
mp.useGPU = false;
mp.flagPlot = false;

%--General
mp.centering = 'pixel';

%--Method of computing core throughput:
% - 'HMI' for energy within half-max isophote divided by energy at telescope pupil
% - 'EE' for encircled energy within a radius (mp.thput_radius) divided by energy at telescope pupil
mp.thput_metric = 'EE'; 
mp.thput_radius = 0.7; %--photometric aperture radius [lambda_c/D]. Used ONLY for 'EE' method.
mp.thput_eval_x = 7; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput

%--Where to shift the source to compute the intensity normalization value.
mp.source_x_offset_norm = 7;  % x location [lambda_c/D] in dark hole at which to compute intensity normalization
mp.source_y_offset_norm = 0;  % y location [lambda_c/D] in dark hole at which to compute intensity normalization

%% Bandwidth and Wavelength Specs

mp.lambda0 = 550e-9;    %--Central wavelength of the whole spectral bandpass [meters]
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

% mp.estimator options:
% - 'perfect' for exact numerical answer from full model
% - 'pairwise' or 'pairwise-square' for pairwise probing in a square region
% centered on the star
% - 'pairwise-rect' for pairwise probing in the specified rectangular regions for
%    one or more stars
% - 'scc' for self-coherent camera
mp.estimator = 'scc';

% %--New variables for pairwise probing estimation:
% mp.est.probe = Probe; % initialize object
% mp.est.probe.Npairs = 3;     % Number of pair-wise probe PAIRS to use.
% mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
% mp.est.probe.radius = 12;    % Max x/y extent of probed region [lambda/D].
% mp.est.probe.xOffset = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
% mp.est.probe.yOffset = 0;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
% mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
% mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.

%% Wavefront Control: General

%--Whether to perform a model-based (instead of empirical) grid search for the controller
mp.ctrl.flagUseModel = true; 

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppress with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)

%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = 2:3; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = [2,3; 3,4; 4,5]; 

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -5:1:-2; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
% % mp.ctrl.dm9regfacVec = 1;        %--Additional regularization factor applied to DM9
   
%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1;
mp.dm2.weight = 1;

%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
mp.controller = 'gridsearchEFC';

% % % GRID SEARCH EFC DEFAULTS     
%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 5; %--Number of estimation+control iterations to perform
mp.dm_ind = [1]; %--Which DMs to use


%% Deformable Mirrors: Influence Functions
%--Influence Function Options:
% - 'influence_dm5v2.fits' for one type of Xinetics DM
% - 'influence_BMC_2kDM_400micron_res10.fits' for BMC 2k DM
% - 'influence_BMC_kiloDM_300micron_res10_spline.fits' for BMC kiloDM
mp.dm1.inf_fn = 'influence_BMC_2kDM_400micron_res10.fits';
mp.dm2.inf_fn = 'influence_BMC_2kDM_400micron_res10.fits';

mp.dm1.dm_spacing = 400e-6; %--User defined actuator pitch [meters]
mp.dm2.dm_spacing = 400e-6; %--User defined actuator pitch [meters]

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';

%% Deformable Mirrors: Optical Layout Parameters

%--DM1 parameters
mp.dm1.Nact = 24;               % # of actuators across DM array
mp.dm1.VtoH = 5e-9*ones(mp.dm1.Nact);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = (mp.dm1.Nact/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = (mp.dm1.Nact/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]
mp.dm1.basisType = 'fourier';  % basis set for control. 'actuator' or 'fourier'

%--DM2 parameters
mp.dm2.Nact = 24;               % # of actuators across DM array
mp.dm2.VtoH = 5e-9*ones(mp.dm1.Nact);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;                % clocking of DM surface [degrees]
mp.dm2.xc = (mp.dm2.Nact/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = (mp.dm2.Nact/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm2.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]
mp.dm2.basisType = 'fourier';  % basis set for control. 'actuator' or 'fourier'

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.dm1.Dstop = mp.dm1.Nact*mp.dm1.dm_spacing;  %--Diameter of iris [meters]
mp.flagDM2stop = false;  %--Whether to apply an iris or not
mp.dm2.Dstop = mp.dm1.Nact*mp.dm1.dm_spacing;   %--Diameter of iris [meters]

%--DM separations
mp.d_P2_dm1 = 0;        % distance (along +z axis) from P2 pupil to DM1 [meters]
mp.d_dm1_dm2 = 0.20;   % distance between DM1 and DM2 [meters]

%% Optical Layout: All models

%--Key Optical Layout Choices
mp.flagSim = true;      %--Simulation or not
mp.layout = 'Fourier';  %--Which optical layout to use
mp.coro = 'vortex';

%--Final Focal Plane Properties
mp.Fend.res = 6; %--Sampling [ pixels per lambda0/D]
% mp.Fend.FOV = 11; %--half-width of the field of view in both dimensions [lambda0/D]
mp.Fend.clockAngDeg = 0; % Clocking angle of the dark hole region

%--Correction and scoring region definition
mp.Fend.corr.Rin = 2;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 5;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 2;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 5;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 180;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'bottom'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'
mp.Fend.shape = 'D';

%% SCC

mp.Fend.Nxi = 150;
mp.Fend.Neta = 150;

% % Jacobian
% mp.path.jac = % Define path to saved out Jacobians. Default if left empty is 'falco-matlab/data/jac/'
% mp.jac.fn = 'jac_scc_test.mat'; % Name of the Jacobian file to save or that is already saved. The path to this file is set by mp.path.jac.
% mp.relinItrVec = 1;  %--Which correction iterations at which to re-compute the control Jacobian. Make empty to reload fn_jac.

mp.scc.modeCoef = 2;  % Gain coeficient to apply to the normalized DM basis set for the empirical SCC Jacobian calculation.
mp.scc.butterworth_exponent = 3;
% mp.scc.pupil_center_row = 46;
% mp.scc.pupil_center_col = 94;
% mp.scc.pupil_subwindow_diameter = 30;

mp.dm1.fourier_spacing = 1; % Center-to-center spacing between Fourier modes in the focal plane. [lambda/D]
mp.dm1.fourier_gridType = 'hex';  % Options: 'hex' or 'square'. 'hex' has a denser packing
xiMin = mp.Fend.corr.Rin-1;
clocking = mp.Fend.clockAngDeg - 90; % -90 for 'bottom' dark hole.
[mp.dm1.fourier_basis_xis , mp.dm1.fourier_basis_etas] = falco_choose_fourier_locations_polar(...
    mp.dm1.Nact/2, mp.dm1.fourier_spacing, mp.dm1.fourier_gridType, xiMin, mp.Fend.corr.Rout+1, mp.Fend.corr.ang, clocking, xiMin);

mp.dm2.fourier_basis_xis = [];
mp.dm2.fourier_basis_etas = [];

mp.dm1.Nactbeam = (mp.dm1.Nact-2); % Number of actuators across the beam (approximate)
mp.dm2.Nactbeam = (mp.dm2.Nact-2); % Number of actuators across the beam (approximate)


%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = (mp.dm1.Nact-2)*mp.dm1.dm_spacing;
mp.P3.D = mp.P2.D;
mp.P4.D = mp.P2.D;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 250;
mp.P2.compact.Nbeam = mp.P1.compact.Nbeam;
mp.P3.compact.Nbeam = mp.P1.compact.Nbeam;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam;  % P4 must be the same as P1 for Vortex. 

%--Number of re-imaging relays between pupil planes in compact model. Needed
%to keep track of 180-degree rotations.
mp.Nrelay1to2 = 0;
mp.Nrelay2to3 = 0;
mp.Nrelay3to4 = 0;
mp.NrelayFend = 0; %--How many times to rotate the final image by 180 degrees

% mp.F3.compact.res = 6; % sampling of FPM for compact model [pixels per lambda0/D]

%% Optical Layout: Full Model 

%--Focal Lengths
% mp.fl = 1; 

%--Pupil Plane Resolutions
mp.P1.full.Nbeam = mp.P1.compact.Nbeam;
mp.P2.full.Nbeam = mp.P1.full.Nbeam;
mp.P3.full.Nbeam = mp.P1.full.Nbeam;
mp.P4.full.Nbeam = mp.P1.full.Nbeam;  % P4 must be the same as P1 for Vortex. 

% mp.F3.full.res = 6;    % sampling of FPM for full model [pixels per lambda0/D]

%% Entrance Pupil (P1) Definition and Generation

mp.P1.D = 4.00; %--circumscribed telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.

% Inputs common to both the compact and full models
clear inputs
inputs.OD = 1.00;

% Full model
inputs.Nbeam = mp.P4.full.Nbeam;
inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
mp.P1.full.mask = falco_gen_pupil_Simple(inputs); 

% Compact model
inputs.Nbeam = mp.P1.compact.Nbeam;
inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
mp.P1.compact.mask = falco_gen_pupil_Simple(inputs); 


%% "Apodizer" (P3) Definition and Generation
mp.flagApod = false;    %--Whether to use an apodizer or not in the FALCO models.


%% Lyot stop (P4) Definition and Generation

clear inputs
DbeamLyot = 9.24e-3; % [meters]
radialOffset = 12.6e-3/DbeamLyot; % [beam diameters]
radiusPinhole = (0.350e-3/2)/DbeamLyot; % [beam diameters]
pinholeRotation = -60; % [degrees]
maxOffset = radialOffset*max([abs(sind(pinholeRotation)), abs(cosd(pinholeRotation))]);
NpadCompact = ceil_even(mp.P4.compact.Nbeam*(2*maxOffset + 2*radiusPinhole) + 1) + 2;
NpadFull = ceil_even(mp.P4.full.Nbeam*(2*maxOffset + 2*radiusPinhole) + 1) + 2;

% Pinhole
inputs.centering = mp.centering;
inputs.OD = 2*radiusPinhole;
inputs.xShear = radialOffset*cosd(pinholeRotation);
inputs.yShear = radialOffset*sind(pinholeRotation);

inputs.Nbeam = mp.P4.compact.Nbeam;
inputs.Npad = NpadCompact;
pinholeCompact = falco_gen_pupil_Simple(inputs); 

inputs.Nbeam = mp.P4.full.Nbeam;
inputs.Npad = NpadFull;
pinholeFull = falco_gen_pupil_Simple(inputs); 

% Main, centered Lyot stop opening
% Inputs common to both the compact and full models
clear inputs
inputs.centering = mp.centering;
inputs.OD = 7.5e-3/DbeamLyot; % IACT values

% Compact model
inputs.Nbeam = mp.P4.compact.Nbeam;
inputs.Npad = NpadCompact;
mp.P4.compact.mask = pinholeCompact + falco_gen_pupil_Simple(inputs); 

% Full model
inputs.Nbeam = mp.P4.full.Nbeam;
inputs.Npad = NpadFull;
mp.P4.full.mask = pinholeFull + falco_gen_pupil_Simple(inputs); 

mp.P4.full.maskWithoutPinhole = falco_gen_pupil_Simple(inputs);
mp.P4.full.maskWithPinhole = pinholeFull + falco_gen_pupil_Simple(inputs); 


figure(1); imagesc(mp.P4.compact.mask); axis xy equal tight; colorbar;
figure(2); imagesc(mp.P4.full.mask); axis xy equal tight; colorbar;


%% Vortex Specific Values %%

mp.F3.VortexCharge = 6; %--Charge of the vortex mask

mp.F3.compact.res = 4; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m
mp.F3.full.res = 4; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m


%% Aberrations

indsZnoll = 30:2:200;
Nz = length(indsZnoll);
ZmapCubeCompact = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam, mp.centering, indsZnoll);
ZmapCubeFull = falco_gen_norm_zernike_maps(mp.P1.full.Nbeam, mp.centering, indsZnoll);

mp.P1.compact.E = exp(1j*10e-9/mp.lambda0*sum(ZmapCubeCompact, 3));
mp.P1.full.E = exp(1j*10e-9/mp.lambda0*sum(ZmapCubeFull, 3));

