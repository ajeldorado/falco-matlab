
% %--Initialize some structures if they don't already exist

%% Misc

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%--Special Computational Settings
mp.flagParfor = true;
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
mp.thput_eval_x = 7; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput

%--Where to shift the source to compute the intensity normalization value.
mp.source_x_offset_norm = 7;  % x location [lambda_c/D] in dark hole at which to compute intensity normalization
mp.source_y_offset_norm = 0;  % y location [lambda_c/D] in dark hole at which to compute intensity normalization

%% Bandwidth and Wavelength Specs

mp.lambda0 = 575e-9;    %--Central wavelength of the whole spectral bandpass [meters]
mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 3;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 3;% 3 or 7  %--Number of wavelengths to used to approximate an image in each sub-bandpass

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

mp.ctrl.flagUseModel = true; %--Whether to perform a model-based (vs empirical) grid search for the controller

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppress with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*[1]; %1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is not used; it is a placeholder for correct indexing of other modes.)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = [2,3]; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = ...
                [3., 4.;
                4.,5.;    
                5., 8.;
                8., 9.];  

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
% % mp.ctrl.dm9regfacVec = 1;        %--Additional regularization factor applied to DM9
   
%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1.;
mp.dm2.weight = 1.;

%--Voltage range restrictions
mp.dm1.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.dm2.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.maxAbsdV = 1000;     %--Max +/- delta voltage step for each actuator for DMs 1 and 2 [volts] %--NOT ENFORCED YET

%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'SM-CVX' for constrained EFC using CVX. --> DEVELOPMENT ONLY

% % % % GRID SEARCH EFC DEFAULTS     
% %--WFSC Iterations and Control Matrix Relinearization
% mp.Nitr = 30; %--Number of estimation+control iterations to perform
% mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
% mp.dm_ind = [1 2]; %--Which DMs to use

% % PLANNED SEARCH EFC DEFAULTS     
mp.dm_ind = [1 2 ]; % vector of DMs used in controller at ANY time (not necessarily all at once or all the time). 
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

% mp.ctrl.sched_mat = [...
%     repmat([1,-2.5,12,1,0],[4,1]);...
%     ];

% mp.ctrl.sched_mat = [...
%     [0,0,0,1,0];...
%     repmat([1,1j,12,0,1],[4,1]);...
%     repmat([1,1j-1,12,0,1],[25,1]);...
%     repmat([1,1j,12,0,1],[1,1]);...
%     ];
% [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

mp.controller = 'plannedEFC';
mp.ctrl.sched_mat = [...
    [0,0,0,1,0];...
    repmat([1,1j,12,0,1],[4,1]);...   %--Optimal beta
    repmat([1,1j-2,12,0,1],[1,1]);... %--Beta kick
    repmat([1,1j,12,0,1],[9,1]);...   %--Optimal beta
    repmat([1,1j-2,12,0,1],[1,1]);... %--Beta kick
    repmat([1,1j,12,0,1],[15,1]);...  %--Optimal beta
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
mp.dm1.ytilt = 5.7;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = 23.5;%(48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = 23.5;%(48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--DM2 parameters
mp.dm2.Nact = 48;               % # of actuators across DM array
mp.dm2.VtoH = 1e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 5.7;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;              % clocking of DM surface [degrees]
mp.dm2.xc = 23.5;%(48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = 23.5;%(48/2 - 1/2);       % y-center location of DM surface [actuator widths]
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
mp.layout = 'wfirst_phaseb_proper';  %--Which optical layout to use. 'wfirst_phaseb_proper' or 'wfirst_phaseb_simple'
mp.coro = 'HLC';
mp.flagApod = false;    %--Whether to use an apodizer or not
mp.flagDMwfe = false;  %--Whether to use BMC DM quilting maps

%--Final Focal Plane Properties
mp.Fend.res = 3;%(730/660)*2.; %--Sampling [ pixels per lambda0/D]
mp.Fend.FOV = 10.; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin = 2.7;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 9.0;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 3;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 9;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 180;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'both'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1.; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = 46.3e-3; %46.2987e-3;
mp.P3.D = 46.3e-3; %46.2987e-3;
mp.P4.D = 46.3e-3; %46.2987e-3;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 309;
mp.P2.compact.Nbeam = 309;
mp.P3.compact.Nbeam = 309;
mp.P4.compact.Nbeam = 309;


%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;
mp.NrelayFend = 1; %--How many times to rotate the final image by 180 degrees

%--FPM resolution
mp.F3.compact.res = 2048./309.;    % sampling of FPM for compact model [pixels per lambda0/D]


%% Optical Layout: Full Model 

mp.full.data_dir = '/Users/ajriggs/Repos/proper-models/wfirst_cgi/data_phaseb/'; % mask design data path
mp.full.cor_type = 'hlc'; %   'hlc', 'spc', or 'none' (none = clear aperture, no coronagraph)

% %--Pupil Plane Resolutions
mp.P1.full.Nbeam = 309;
mp.P1.full.Narr = 310;

mp.full.output_dim = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D
mp.full.flagPROPER = true; %--Whether the full model is a PROPER prescription

mp.full.pol_conds = 10;% [-2,-1,1,2]; %--Which polarization states to use when creating an image.
mp.full.polaxis = 10;                %   polarization condition (only used with input_field_rootname)
mp.full.use_errors = true;

mp.full.zindex = 4;
mp.full.zval_m = 0.19e-9;
mp.full.use_hlc_dm_patterns = false; % whether to use design WFE maps for HLC.
mp.full.lambda0_m = mp.lambda0;
mp.full.input_field_rootname = '';	%   rootname of files containing aberrated pupil
mp.full.use_dm1 = 0;                %   use DM1? 1 or 0
mp.full.use_dm2 = 0;                %   use DM2? 1 or 0
mp.full.dm_sampling_m = 0.9906e-3;     %   actuator spacing in meters; default is 1 mm
mp.full.dm1_xc_act = 23.5;          %   for 48x48 DM, wavefront centered at actuator intersections: (0,0) = 1st actuator center
mp.full.dm1_yc_act = 23.5;
mp.full.dm1_xtilt_deg = 0;   		%   tilt around X axis
mp.full.dm1_ytilt_deg = 5.7;		%   effective DM tilt in deg including 9.65 deg actual tilt and pupil ellipticity
mp.full.dm1_ztilt_deg = 0;
mp.full.dm2_xc_act = 23.5;		
mp.full.dm2_yc_act = 23.5;
mp.full.dm2_xtilt_deg = 0;   
mp.full.dm2_ytilt_deg = 5.7;
mp.full.dm2_ztilt_deg = 0;
mp.full.use_fpm  = 1;
mp.full.fpm_axis = 'p';             %   HLC FPM axis: '', 's', 'p'

mp.full.dm1.flatmap = 0;%0.5*fitsread('errors_polaxis10_dm.fits');
mp.full.dm2.flatmap = 0;%0.5*fitsread('errors_polaxis10_dm.fits');


%% Mask Definitions

%--Pupil definition
mp.whichPupil = 'WFIRST180718';
mp.P1.IDnorm = 0.303; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.D = 2.3631; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
mp.P1.Dfac = 1; %--Factor scaling inscribed OD to circumscribed OD for the telescope pupil.

% %--Lyot stop shape
% mp.P4.IDnorm = 0.50; %--Lyot stop ID [Dtelescope]
% mp.P4.ODnorm = 0.80; %--Lyot stop OD [Dtelescope]
% % mp.P4.ang = 90;      %--Lyot stop opening angle [degrees]
% mp.P4.wStrut = 0.036;    %--Lyot stop strut width [pupil diameters]

% %--FPM size
% mp.F3.Rin = 2.8;   % inner hard-edge radius of the focal plane mask [lambda0/D]. Needs to be <= mp.F3.Rin 
% mp.F3.Rout = Inf;   % radius of outer opaque edge of FPM [lambda0/D]
% mp.F3.ang = 180;    % on each side, opening angle [degrees]

%--Whether to generate or load various masks: compact model
mp.compact.flagGenPupil = false;  
mp.compact.flagGenFPM = false;  
mp.compact.flagGenLS = false;

mp.P1.compact.mask = fitsread([mp.full.data_dir 'hlc_20190210/' 'run461_pupil.fits']);
mp.P1.compact.mask = padOrCropEven(mp.P1.compact.mask,310);

mp.P4.compact.mask = fitsread([mp.full.data_dir 'hlc_20190210/' 'run461_lyot.fits']);
mp.P4.compact.mask = padOrCropEven(mp.P4.compact.mask,310);

% mp.dm1.wfe = fitsread([mp.full.data_dir 'hlc_20190210/' 'run461_dm1wfe.fits']);
% mp.dm2.wfe = fitsread([mp.full.data_dir 'hlc_20190210/' 'run461_dm2wfe.fits']);
