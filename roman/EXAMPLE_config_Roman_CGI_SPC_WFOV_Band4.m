% %--Initialize some structures if they don't already exist

%% Path to data needed by PROPER model

mp.full.data_dir = '/Users/ajriggs/Documents/Sim/cgi/public/roman_phasec_v1.2.4/phasec_data/';


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

%--Method of computing core throughput:
% - 'HMI' for energy within half-max isophote divided by energy at telescope pupil
% - 'EE' for encircled energy within a radius (mp.thput_radius) divided by energy at telescope pupil
mp.thput_metric = 'HMI'; 
mp.thput_radius = 0.7; %--photometric aperture radius [lambda_c/D]. Used ONLY for 'EE' method.
mp.thput_eval_x = 13; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput

%--Where to shift the source to compute the intensity normalization value.
mp.source_x_offset_norm = 13;  % x location [lambda_c/D] in dark hole at which to compute intensity normalization
mp.source_y_offset_norm = 0;  % y location [lambda_c/D] in dark hole at which to compute intensity normalization

%% Bandwidth and Wavelength Specs

mp.lambda0 = 825e-9;   %--Central wavelength of the whole spectral bandpass [meters].
mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 3;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 3;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing in the specified rectangular regions for
%    one or more stars
% - 'pwp-bp-square' for pairwise probing with batch process estimation in a
% square region for one star [original functionality of 'pwp-bp' prior to January 2021]
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT TESTED YET]
mp.estimator = 'perfect';

%--New variables for pairwise probing estimation:
mp.est.probe.Npairs = 3;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 21;    % Max x/y extent of probed region [lambda/D].
mp.est.probe.xOffset = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.yOffset = 14;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.

%% Wavefront Control: General

mp.ctrl.flagUseModel = true; %--Whether to perform a model-based (vs empirical) grid search for the controller

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppress with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*[1]; %1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is not used; it is a placeholder for correct indexing of other modes.)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = [];%2:3;%2:6; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = [];%...
%                 [3., 4.;...
%                 4., 8.;
%                 8., 9.];  

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
% % mp.ctrl.dm9regfacVec = 1;        %--Additional regularization factor applied to DM9
   
%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1.;
mp.dm2.weight = 1.;

%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
mp.controller = 'gridsearchEFC';

% % % GRID SEARCH EFC DEFAULTS     
%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 5; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.dm_ind = [1 2]; %--Which DMs to use

% % % PLANNED SEARCH EFC DEFAULTS     
% mp.dm_ind = [1 2 ]; % vector of DMs used in controller at ANY time (not necessarily all at once or all the time). 
% mp.ctrl.dmfacVec = 1;
% %--CONTROL SCHEDULE. Columns of mp.ctrl.sched_mat are: 
%     % Column 1: # of iterations, 
%     % Column 2: log10(regularization), 
%     % Column 3: which DMs to use (12, 128, 129, or 1289) for control
%     % Column 4: flag (0 = false, 1 = true), whether to re-linearize
%     %   at that iteration.
%     % Column 5: flag (0 = false, 1 = true), whether to perform an
%     %   EFC parameter grid search to find the set giving the best
%     %   contrast .
%     % The imaginary part of the log10(regularization) in column 2 is
%     %  replaced for that iteration with the optimal log10(regularization)
%     % A row starting with [0, 0, 0, 1...] is for relinearizing only at that time
% 
% mp.ctrl.sched_mat = [...
%     [0,0,0,1,0];
%     repmat([1,1j,12,0,1],[10,1]);...
%     ];
% 
% % mp.ctrl.sched_mat = [...
% %     repmat([1,1j,12,1,1],[4,1]);...
% %     repmat([1,1j-1,12,1,1],[25,1]);...
% %     repmat([1,1j,12,1,1],[1,1]);...
% %     ];
% [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);


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
mp.dm1.orientation = 'rot180'; % Change to mp.dm1.V orientation before generating DM surface. Options: rot0, rot90, rot180, rot270, flipxrot0, flipxrot90, flipxrot180, flipxrot270
mp.dm1.Nact = 48;               % # of actuators across DM array
mp.dm1.VtoH = 4e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 9.65;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = 23.5;%(48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = 23.5;%(48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]
mp.dm1.Vmin = 0; %--Min allowed absolute voltage command
mp.dm1.Vmax = 100; %--Max allowed absolute voltage command
mp.dm1.dVnbrLat = 50; % max voltage difference allowed between laterally-adjacent DM actuators
mp.dm1.dVnbrDiag = 75; % max voltage difference allowed between diagonally-adjacent DM actuators
mp.dm1.facesheetFlatmap = 50 * ones(mp.dm1.Nact, mp.dm1.Nact); %--Voltage map that produces a flat DM2 surface. Used when enforcing the neighbor rule.

%--DM2 parameters
mp.dm2.orientation = 'rot180'; % Change to mp.dm2.V orientation before generating DM surface. Options: rot0, rot90, rot180, rot270, flipxrot0, flipxrot90, flipxrot180, flipxrot270
mp.dm2.Nact = 48;               % # of actuators across DM array
mp.dm2.VtoH = 4e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 9.65;%               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;              % clocking of DM surface [degrees]
mp.dm2.xc = 23.5;%(48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = 23.5;%(48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm2.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]
mp.dm2.Vmin = 0; %--Min allowed absolute voltage command
mp.dm2.Vmax = 100; %--Max allowed absolute voltage command
mp.dm2.dVnbrLat = 50; % max voltage difference allowed between laterally-adjacent DM actuators
mp.dm2.dVnbrDiag = 75; % max voltage difference allowed between diagonally-adjacent DM actuators
mp.dm2.facesheetFlatmap = 50 * ones(mp.dm2.Nact, mp.dm2.Nact); %--Voltage map that produces a flat DM2 surface. Used when enforcing the neighbor rule.

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.dm1.Dstop = 100e-3;  %--Diameter of iris [meters]
mp.flagDM2stop = true;  %--Whether to apply an iris or not
mp.dm2.Dstop = 51.5596e-3;   %--Diameter of iris [meters]

%--DM separations
mp.d_P2_dm1 = 0;        % distance (along +z axis) from P2 pupil to DM1 [meters]
mp.d_dm1_dm2 = 1.000;   % distance between DM1 and DM2 [meters]


%% Optical Layout: All models

%--Key Optical Layout Choices
mp.flagSim = true;      %--Simulation or not
mp.layout = 'roman_phasec_proper';  %--Which optical layout to use
mp.coro = 'SPLC';
mp.flagRotation = false;
mp.flagApod = true;    %--Whether to use an apodizer or not
mp.flagDMwfe = false;  %--Whether to apply DM aberration maps in FALCO models

%--Final Focal Plane Properties
mp.Fend.res = 3.30; %--Sampling [ pixels per lambda0/D]. 825/500*2
mp.Fend.FOV = 22.0; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin = 5.6;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 20.4;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 6.0;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 20.0;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 180;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'lr'; %--Which side(s) for correction: 'left', 'right', 'top', 'up', 'bottom', 'down', 'lr', 'rl', 'leftright', 'rightleft', 'tb', 'bt', 'ud', 'du', 'topbottom', 'bottomtop', 'updown', 'downup'
mp.Fend.clockAngDeg = 0; %--Amount to rotate the dark hole location


%% Optical Layout: Full Model 

mp.full.cor_type = 'spc-wide_band4';

mp.full.flagPROPER = true; %--Whether the full model is a PROPER prescription

% %--Pupil Plane Resolutions
mp.P1.full.Nbeam = 1000;
mp.P1.full.Narr = 1002;

mp.full.output_dim = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D

mp.full.pol_conds = 10;% [-2,-1,1,2]; %--Which polarization states to use when creating an image.
mp.full.polaxis = 10;                %   polarization condition (only used with input_field_rootname)
mp.full.use_errors = true;

mp.full.dm1.flatmap = fitsread('spc_wide_band4_flattened_dm1.fits');
mp.full.dm2.flatmap = fitsread('spc_wide_band4_flattened_dm2.fits');

mp.dm1.biasMap = 50 + mp.full.dm1.flatmap./mp.dm1.VtoH; %--Bias voltage. Needed prior to WFSC to allow + and - voltages. Total voltage is mp.dm1.biasMap + mp.dm1.V
mp.dm2.biasMap = 50 + mp.full.dm2.flatmap./mp.dm2.VtoH; %--Bias voltage. Needed prior to WFSC to allow + and - voltages. Total voltage is mp.dm2.biasMap + mp.dm2.V

%% Optical Layout: Compact Model (and Jacobian Model)

%--Focal Lengths
mp.fl = 1.; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = 46.3e-3;
mp.P3.D = 46.3e-3;
mp.P4.D = 46.3e-3;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 300;
mp.P2.compact.Nbeam = 300;
mp.P3.compact.Nbeam = 300;
mp.P4.compact.Nbeam = 120;

%--Shaped Pupil Mask: Load and downsample.
% mp.SPname = 'SPC-20190130';
SP0 = fitsread([mp.full.data_dir filesep 'spc_20200610_wfov' filesep 'SPM_SPC-20200610_1000_rounded9_gray.fits']);
SP0 = pad_crop(SP0, 1001);
SP0 = rot90(SP0, 2);
% % SP0(2:end,2:end) = rot90(SP0(2:end,2:end),2);

SP1 = falco_filtered_downsample(SP0, mp.P3.compact.Nbeam/mp.P1.full.Nbeam, mp.centering);
mp.P3.compact.mask = pad_crop(SP1, ceil_even(max(size(SP1))));

%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
% NOTE: All these relays are ignored if mp.flagRotation == false.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;
mp.NrelayFend = 1; %--How many times to rotate the final image by 180 degrees


%% Mask Definitions

%--Pupil definition
mp.P1.IDnorm = 0.303; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.D = 2.3631; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
mp.P1.Dfac = 1; %--Factor scaling inscribed OD to circumscribed OD for the telescope pupil.
changes.flagRot180 = true;
mp.P1.compact.mask = falco_gen_pupil_Roman_CGI_20200513(mp.P1.compact.Nbeam, mp.centering, changes);

%--Lyot stop 
mp.P4.IDnorm = 0.36; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.91; %--Lyot stop OD [Dtelescope]
wStrut = 3.2/100; % Lyot stop strut width [pupil diameters]
rocFilletLS = 0.02; % [pupil diameters]
upsampleFactor = 100; %--Lyot anti-aliasing value
mp.P4.compact.mask = falco_gen_Roman_CGI_lyot_stop_symm_fillet(mp.P4.compact.Nbeam, mp.P4.IDnorm, mp.P4.ODnorm, wStrut, rocFilletLS, upsampleFactor, mp.centering);

% FPM parameters
mp.F3.compact.res = 3;
inputs.pixresFPM = mp.F3.compact.res; % [pixels per lambda0/D]
inputs.rhoInner = 5.6; % [lambda0/D]
inputs.rhoOuter = 20.4; % [lambda0/D]
mp.F3.compact.mask = falco_gen_annular_FPM(inputs);
