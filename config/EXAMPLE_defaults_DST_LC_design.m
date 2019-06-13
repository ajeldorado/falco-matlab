
% %--Initialize some structures if they don't already exist

%% Misc

%--Record Keeping
mp.SeriesNum = 867;
mp.TrialNum = 5309;

%--Special Computational Settings
mp.flagParfor = false;
mp.useGPU = false;
mp.flagPlot = false;

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

mp.lambda0 = 550e-9;    %--Central wavelength of the whole spectral bandpass [meters]
mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 5;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing with batch process estimation
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT TESTED YET]
% - 'pwp-iekf' for pairwise probing with iterated extended Kalman filter  [NOT AVAILABLE YET]
mp.estimator = 'pwp-bp';

%--New variables for pairwise probing estimation:
mp.est.probe.Npairs = 3;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 12;    % Max x/y extent of probed region [actuators].
mp.est.probe.offsetX = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.offsetY = 10;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
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
mp.eval.indsZnoll = 2:3; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = [3, 4;...
              4, 8];  

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
% % mp.ctrl.dm9regfacVec = 1;        %--Additional regularization factor applied to DM9
   
%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1;
mp.dm2.weight = 1;

%--Voltage range restrictions
mp.dm1.maxAbsV = 1000;
mp.dm2.maxAbsV = 1000;
mp.maxAbsdV = 1000;  %--Max +/- delta voltage step for each actuator for DMs 1 and 2


%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'SM-CVX' for constrained EFC using CVX. --> DEVELOPMENT ONLY
mp.controller = 'plannedEFC';

% % % % GRID SEARCH EFC DEFAULTS     
% %--WFSC Iterations and Control Matrix Relinearization
% mp.Nitr = 20; %--Number of estimation+control iterations to perform
% mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
% 
% mp.dm_ind = [1 2]; %--Which DMs to use


% % % PLANNED SEARCH EFC DEFAULTS     
mp.dm_ind = [1 2]; % vector of DMs used in controller at ANY time (not necessarily all at once or all the time). 
mp.ctrl.dmfacVec = 1;
%--CONTROL SCHEDULE. Columns of mp.ctrl.sched_mat are: 
    % Column 1: # of iterations, 
    % Column 2: log10(regularization), 
    % Column 3: which DMs to use (12, 128, 129, or 1289) for control
    % Column 4: flag (0 = false, 1 = true), whether to re-linearize
    %   at that iteration.
    % Column 5: flag (0 = false, 1 = true), whether to perform an
    %   EFC parameter grid search to find the set giving the best
    %   contrast .
    % The imaginary part of the log10(regularization) in column 2 is
    %  replaced for that iteration with the optimal log10(regularization)
    % A row starting with [0, 0, 0, 1...] is for relinearizing only at that time
mp.ctrl.sched_mat = [...
    repmat([1,1j,12,1,1],[4,1]);...
    repmat([1,1j-1,12,1,1],[25,1]);...
    repmat([1,1j,12,1,1],[1,1]);...
    ];
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%% Deformable Mirrors: Influence Functions
%--Influence Function Options:
% - 'influence_dm5v2.fits' for one type of Xinetics DM
% - 'influence_BMC_2kDM_400micron_res10.fits' for BMC 2k DM
% - 'influence_BMC_kiloDM_300micron_res10_spline.fits' for BMC kiloDM

mp.dm1.inf_fn = 'influence_dm5v2.fits';
mp.dm2.inf_fn = 'influence_dm5v2.fits';

mp.dm1.dm_spacing = 1e-3; %--User defined actuator pitch
mp.dm2.dm_spacing = 1e-3; %--User defined actuator pitch

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';

%% Deformable Mirrors: Optical Layout Parameters

%--DM1 parameters
mp.dm1.Nact = 48;               % # of actuators across DM array
mp.dm1.VtoH = 1*1e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = (48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = (48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--DM2 parameters
mp.dm2.Nact = 48;               % # of actuators across DM array
mp.dm2.VtoH = 1*1e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;                % clocking of DM surface [degrees]
mp.dm2.xc = (48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = (48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm2.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.dm1.Dstop = 100e-3;  %--Diameter of iris [meters]
mp.flagDM2stop = true;  %--Whether to apply an iris or not
mp.dm2.Dstop = 50e-3;   %--Diameter of iris [meters]

%--DM separations
mp.d_P2_dm1 = 0;        % distance (along +z axis) from P2 pupil to DM1 [meters]
mp.d_dm1_dm2 = 1.000;   % distance between DM1 and DM2 [meters]


%% Optical Layout: All models

%--Key Optical Layout Choices
mp.flagSim = true;      %--Simulation or not
mp.layout = 'Fourier';  %--Which optical layout to use
mp.coro = 'LC';
mp.flagApod = false;    %--Whether to use an apodizer or not

%--Final Focal Plane Properties
mp.Fend.res = 3; %--Sampling [ pixels per lambda0/D]
mp.Fend.FOV = 11; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin = 2.8;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 2.8;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 10;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 180;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'both'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = 46.3e-3;
mp.P3.D = 46.3e-3;
mp.P4.D = 46.3e-3;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 300;
mp.P2.compact.Nbeam = 300;
mp.P3.compact.Nbeam = 300;
mp.P4.compact.Nbeam = 300;

%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;
mp.NrelayFend = 0; %--How many times to rotate the final image by 180 degrees

mp.F3.compact.res = 6; % sampling of FPM for compact model [pixels per lambda0/D]

%% Optical Layout: Full Model 

%--Focal Lengths
% mp.fl = 1; 

%--Pupil Plane Resolutions
mp.P1.full.Nbeam = 300;
mp.P2.full.Nbeam = 300;
mp.P3.full.Nbeam = 300;
mp.P4.full.Nbeam = 300;

mp.F3.full.res = 6;    % sampling of FPM for full model [pixels per lambda0/D]

%% Mask Definitions

%--Pupil definition
mp.whichPupil = 'Simple';
mp.P1.IDnorm = 0.00; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.ODnorm = 1.00;% Outer diameter of the telescope [diameter]
mp.P1.stretch = 1.00; % factor that stretches the horizontal axis to create elliptical beam 
mp.P1.D = 4; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
mp.P1.Dfac = 1; %--Factor scaling inscribed OD to circumscribed OD for the telescope pupil.
mp.P1.Nstrut = 0;% Number of struts 
mp.P1.angStrut = [];%Array of angles of the radial struts (deg)
mp.P1.wStrut = []; % Width of the struts (fraction of pupil diam.)

%--Lyot stop padding
mp.P4.IDnorm = 47.36/227.86; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 156.21/227.86; %--Lyot stop OD [Dtelescope]
mp.P4.Nstrut = 3;% Number of struts 
mp.P4.angStrut = [90 210 330];%Array of angles of the radial struts (deg)
mp.P4.wStrut = 0.005; % Width of the struts (fraction of pupil diam.)

%--FPM size
mp.F3.Rin = 2.8;    % radius of inner hard edge of the focal plane mask [lambda0/D]
mp.F3.Rout = Inf;   % radius of outer opaque edge of FPM [lambda0/D]
mp.F3.ang = 180;    % on each side, opening angle [degrees]
mp.FPMampFac = 10^(-3.7); % amplitude transmission of the FPM


%% LC-Specific Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


