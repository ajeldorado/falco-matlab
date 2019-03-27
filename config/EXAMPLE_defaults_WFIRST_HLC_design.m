
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
mp.thput_metric = 'HMI'; 
mp.thput_radius = 0.7; %--photometric aperture radius [lambda_c/D]. Used ONLY for 'EE' method.
mp.thput_eval_x = 6; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput

%% Bandwidth and Wavelength Specs

mp.lambda0 = 575e-9;    %--Central wavelength of the whole spectral bandpass [meters]
mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 6;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing with batch process estimation
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT TESTED YET]
% - 'pwp-iekf' for pairwise probing with iterated extended Kalman filter  [NOT AVAILABLE YET]
mp.estimator = 'perfect';

%--New variables for pairwise probing estimation:
mp.est.probe.Npairs = 3;%2;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 12;%20;    % Max x/y extent of probed region [actuators].
mp.est.probe.offsetX = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.offsetY = 14;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.

%--New variables for pairwise probing with a Kalman filter
%  mp.est.ItrStartKF =  %Which correction iteration to start recursive estimate
%  mp.est.tExp =
%  mp.est.num_im =
%  mp.readNoiseStd =
%  mp.peakCountsPerPixPerSec =
%  mp.est.Qcoef =
%  mp.est.Rcoef =

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

%--Voltage range restrictions: general
mp.dm1.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.dm2.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.maxAbsdV = 1000;     %--Max +/- delta voltage step for each actuator for DMs 1 and 2 [volts] %--NOT ENFORCED YET

%--Voltage range restrictions: neighboring actuators
mp.dm1.flagNbrRule = true;
mp.dm1.dVnbr = 150; %--absolute value of max delta voltage between neighbors [volts]
mp.dm2.flagNbrRule = true;
mp.dm2.dVnbr = 150; %--absolute value of max delta voltage between neighbors [volts]

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
% mp.dm_ind = [1 2]; %--Which DMs to use

% % PLANNED SEARCH EFC DEFAULTS     
mp.dm_ind = [1 2 9]; % vector of DMs used in controller at ANY time (not necessarily all at once or all the time). 
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

% SetA = ... %--DMs 1 & 2 for x iterations. Relinearize every iteration.
%     repmat([1, 1j, 12, 1, 1], [20, 1]); 
% SetB = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
%     [0, 0, 0, 1, 0;...
%     10, -3, 129, 0, 0;...
%     5,  -4, 129, 0, 0;...
%     10, -2, 129, 0, 0;...
%     ];
% SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
%     [0, 0, 0, 1, 0;...
%     10, -4, 129, 0, 0;...
%     5,  -5, 129, 0, 0;...
%     10, -2, 129, 0, 0;...
%     ];
% SetD = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
%    [0, 0, 0, 1, 0;...
%    10, -5, 129, 0, 0;...
%    5,  -6, 129, 0, 0;...
%    10, -2, 129, 0, 0;...
%    ];
% SetA2 = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
% SetB2 = [1, -5, 12, 1, 0];
% SetC2 = [1, 1j, 12, 1, 1];
% 
% mp.ctrl.sched_mat = [...
%    repmat(SetA2,[5,1]);...
%    repmat(SetB2,[3,1]);...
%    repmat(SetC2,[3,1]);...
%    ...repmat(SetB,[2,1]);...
%    repmat(SetC,[4,1]);...
%    repmat(SetD,[4,1]);...
%    ];

SetJ = [...
    repmat([1,-5,129,1,1],[4,1]);...
    repmat([1,1j-1,129,1,1],[4,1]);...
    repmat([1,1j,129,1,1],[2,1]);...
    ];
mp.ctrl.sched_mat = [...
    repmat([1,1j,  12,1,1],[5,1]);...
    repmat([1,1j-1,12,1,1],[6,1]);...
    repmat([1,1j,  12,1,1],[1,1]);...
    repmat(SetJ,[10,1]);...
    ];

[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%% Deformable Mirrors: Influence Functions
%--Influence Function Options:
% - 'influence_dm5v2.fits' for one type of Xinetics DM
% - 'influence_BMC_2kDM_400micron_res10.fits' for BMC 2k DM
% - 'influence_BMC_kiloDM_300micron_res10_spline.fits' for BMC kiloDM

mp.dm1.inf_fn = 'influence_dm5v2.fits';
mp.dm2.inf_fn = 'influence_dm5v2.fits';

mp.dm1.dm_spacing = 0.9906e-3;%1e-3; %--User defined actuator pitch
mp.dm2.dm_spacing = 0.9906e-3;%1e-3; %--User defined actuator pitch

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';

%% Deformable Mirrors: Optical Layout Parameters

%--DM1 parameters
mp.dm1.Nact = 48;               % # of actuators across DM array
mp.dm1.VtoH = 1e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 5.83;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = (48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = (48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--DM2 parameters
mp.dm2.Nact = 48;               % # of actuators across DM array
mp.dm2.VtoH = 1e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 5.55;               % for foreshortening. angle of rotation about y-axis [degrees]
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
mp.coro = 'HLC';
mp.flagApod = false;    %--Whether to use an apodizer or not

%--Final Focal Plane Properties
mp.Fend.res = 3; %--Sampling [ pixels per lambda0/D]
mp.Fend.FOV = 12; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin = 2.7;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 2.7;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 10;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 180;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'both'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = 46.2987e-3;
mp.P3.D = 46.2987e-3;
mp.P4.D = 46.2987e-3;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 250;
mp.P2.compact.Nbeam = 250;
mp.P3.compact.Nbeam = 250;
mp.P4.compact.Nbeam = 250;

%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;

%% Optical Layout: Full Model 

%--Focal Lengths
% mp.fl = 1; 


%--Pupil Plane Resolutions
mp.P1.full.Nbeam = 250;
mp.P2.full.Nbeam = 250;
mp.P3.full.Nbeam = 250;
mp.P4.full.Nbeam = 250;

%% Mask Definitions

%--Pupil definition
mp.whichPupil = 'WFIRST180718';
mp.P1.IDnorm = 0.303; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.D = 2.3631; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
mp.P1.Dfac = 1; %--Factor scaling inscribed OD to circumscribed OD for the telescope pupil.

%--Lyot stop padding
mp.P4.wStrut = 3.6/100.; % nominal pupil's value is 76mm = 3.216%
mp.P4.IDnorm = 0.45; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.78; %--Lyot stop OD [Dtelescope]

%--FPM size
mp.F3.Rin = 2.7;    % maximum radius of inner part of the focal plane mask [lambda0/D]
mp.F3.RinA = 2.7;   % inner hard-edge radius of the focal plane mask [lambda0/D]. Needs to be <= mp.F3.Rin 
mp.F3.Rout = Inf;   % radius of outer opaque edge of FPM [lambda0/D]
mp.F3.ang = 180;    % on each side, opening angle [degrees]




%% HLC-Specific Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% FPM Material Properties
mp.aoi = 6.69;%10.04; % Angle of incidence at FPM [deg]
mp.t_Ti_nm = 3.0; %--Static base layer of titanium beneath any nickel [nm]

mp.dt_metal_nm = 1;%0.25;%1/10; %--thickness step size for FPM metal layer (nm)
mp.t_metal_nm_vec = 0:mp.dt_metal_nm:120; %150; %--nickel thickness range and sampling (nm)
mp.dt_diel_nm = 1;%2/10; %--thickness step size for FPM dielectric layer  (nm)
mp.t_diel_nm_vec = 0:mp.dt_diel_nm:900; %--PMGI thickness range and sampling (nm)

%--Number of waves offset from substrate for reference plane. MUST BE MORE THAN MAX THICKNESS OF THE FPM.
mp.FPM.d0fac = 4;

%% DM8: FPM Metal Thickness

mp.dm8.V0coef = 100; % Nominal Nickel layer thickness [nm]

%--DM8 parameters
mp.dm8.Vmin = 0;
mp.dm8.Vmax = 300;


%% DM9: FPM Dielectric thickness

%--DM9 weights and sensitivities: Used by the controller
mp.dm9.weight = 1; % Jacobian weight for the FPM dielectric. Smaller weight makes stroke larger by the inverse of this factor.
mp.dm9.act_sens = 10; %--Change in oomph (E-field sensitivity) of DM9 actuators. Chosen empirically based on how much DM9 actuates during a control step.
mp.dm9.stepFac = 10;%200; %--Adjust the step size in the Jacobian, then divide back out. Used for helping counteract effect of discretization.

%--Starting dielectric thicknesses
mp.t_diel_bias_nm = 0; %--Thickness of starting uniform bias layer of PMGI [nm]. % (Requires an outer stop in reality if >0, but will run without it to see if it gives essentially the same result as the EHLC but faster)
mp.dm9.V0coef = 390; % Nominal PMGI layer thickness [nm] 

%--DM9 influence function options:
% - '3x3'
% - 'Xinetics'
% - '3foldZern'
%--DM9 parameters for 3x3 influence function
mp.dm9.actres = 7; % number of "actuators" per lambda0/D in the FPM's focal plane. On a square actuator array.
mp.dm9.FPMbuffer = -0.5; %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)
mp.dm9.inf0name = '3x3';   % This gives inf0 = 1/4*[1, 2, 1; 2, 4, 2; 1, 2, 1];  

% %%--DM9 parameters for Lanczos3 influence function
% mp.dm9.actres = 8;% % number of "actuators" per lambda0/D in the FPM's focal plane. On a square actuator array.
% mp.dm9.FPMbuffer = 0.2; %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)
% mp.dm9.inf0name = 'Lanczos3'; 
% %--FPM resolution (pixels per lambda0/D) in the compact and full models.
% mp.F3.compact.res = 30;
% mp.F3.full.res = 30;

% %%--DM9 parameters for Xinetics influence function
% mp.dm9.actres = 8;% % number of "actuators" per lambda0/D in the FPM's focal plane. On a square actuator array.
% mp.dm9.FPMbuffer = 0.2; %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)
% mp.dm9.inf0name = 'Xinetics'; 
% %--FPM resolution (pixels per lambda0/D) in the compact and full models.
% mp.F3.compact.res = 30;
% mp.F3.full.res = 30;

%--DM9 parameters
mp.dm9.Vmin = 0;  % minimum thickness of FPM dielectric layer (nm)
mp.dm9.Vmax = 400+mp.dm9.V0coef; % maximum thickness (from one actuator, not of the facesheet) of FPM dielectric layer (nm)

