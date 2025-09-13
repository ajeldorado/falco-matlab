%% DST configuration 

%% Add this file to the list of scripts to copy over
mp.dummy = 1;
nameOfThisFile = mfilename('fullpath'); %[mfilename('fullpath') '.m'];
if ~isfield(mp, 'filesToCopy')
    mp.filesToCopy{1} = nameOfThisFile; % create cell array
else % append to cell array
    Nfn = length(mp.filesToCopy);
    mp.filesToCopy{Nfn+1} = nameOfThisFile;
end

%%

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
% - 'HMI' for energy within half-max isophote divided by energmp.Fend.sidesy at telescope pupil
% - 'EE' for encircled energy within a radius (mp.thput_radius) divided by energy at telescope pupil
mp.thput_metric = 'EE'; 
mp.thput_radius = 0.7; %--photometric aperture radius [lambda_c/D]. Used ONLY for 'EE' method.
mp.thput_eval_x = 6; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput
mp.source_x_offset_norm = 6;
mp.source_y_offset_norm = 0;

%% Bandwidth and Wavelength Specs

mp.lambda0 = 660e-9;    %--Central wavelength of the whole spemp.Fend.sidesctral bandpass [meters]
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing with batch process estimation
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT TESTED YET]
% - 'pwp-iekf' for pairwise probing with iterated extended Kalman filter  [NOT AVAILABLE YET]
mp.estimator = 'pwp-bp-square';

%--New variables for pairwise probing estimation:
mp.est.probe = Probe; 
mp.est.probe.Npairs = 2;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 20;    % Max x/y extent of probed region [actuators].
mp.est.probe.xOffset = 7;    % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.yOffset = -7;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'y';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate] %y for left, x for top
mp.est.probe.gainFudge = 1; % empirical fudge factor to make average probe amplitude match desired value.
mp.est.Ithreshold = 1e-2;

%--New variables for pairwise probing with a Kalman filter
%  mp.est.ItrStartKF =  %Which correction iteration to start recursive estimate
%  mp.est.tExp =
%  mp.est.num_im =
%  mp.readNoiseStd =
%  mp.peakCountsPerPixPerSec =
%  mp.est.Qcoef =
%  mp.est.Rcoef =mp.Fend.sides

%% Wavefront Control: General

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppres with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = []; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = [];  


%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1;
mp.dm2.weight = 1;

%--Voltage range restrictions
mp.maxAbsdV = 50;  %--Max +/- delta voltage step for each actuator for DMs 1 and 2


%% Wavefront Control: Controller Specific
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
mp.dm_ind = 1; %--Which DMs to use
mp.ctrl.log10regVec = -3:-1; %--log10 of the regularization exponents (often called Beta values)

% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schmp.Fend.sidesedule
%  - 'SM-CVX' for constrained EFC using CVX. --> DEVELOPMENT ONLY
% mp.controller = 'gridsearchEFC';
mp.controller = 'gridsearchEFC';

if(strcmp(mp.controller,'gridsearchEFC'))
    % % % GRID SEARCH EFC DEFAULTS     
    %--WFSC Iterations and Control Matrix Relinearization
    mp.Nitr = 10; %--Number of estimation+control iterations to perform
    %mp.relinItrVec = [1,10:10:(mp.Nitr-10)];%1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
    mp.relinItrVec = 1:mp.Nitr; 
    %mp.relinItrVec = 1; 
    %--Grid- or Line-Search Settings
    
elseif(strcmp(mp.controller,'plannedEFC'))

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

end
%% Deformable Mirrors

falco_dst2_dmparams_BMC50E_T2K_20250716; 

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.flagDM2stop = false;  %--Whether to apply an iris or not

% mp.dm2.Dstop = 

%--DM separations
mp.d_P2_dm1 = 0;        % distance (along +z axis) from P2 pupil to DM1 [meters]
mp.d_dm1_dm2 = 0;   % distance between DM1 and DM2 [meters]

%% Optical Layout: All models

%--Key Optical Layout Choices
mp.flagSim = true;      %--Simulation or not
mp.layout = 'Fourier';  %--Which optical layout to use% 

mp.flagApod = false;    %--Whether to use an apodizer or not


%--Final Focal Plane Properties
% mp.Fend.res = 4;%tb.info.pixelPerLamOverD; %--Sampling [ pixels per lambda0/D]
% mp.Fend.FOV = 20;%(tb.sciCam.subwindowsize-2)/(2*mp.Fend.res); %--half-width of the field of view in both dimensions [lambda0/D]
% FALCO adds two pixels back in 

%--Correction and scoring region definition
mp.Fend.score.Rin = 6;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 9; % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 30; % angular opening of dark hole scoring region [degrees]

mp.Fend.corr.Rin = mp.Fend.score.Rin;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout = mp.Fend.score.Rout;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = mp.Fend.score.ang;  % angular opening of dark hole correction region [degrees]

mp.Fend.sides = 'top'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.
% 
%--Pupil Plane Diameters
mp.P2.D = mp.dm1.Nactbeam*mp.dm1.dm_spacing;
mp.P3.D = mp.P2.D;
mp.P4.D = mp.P2.D;
% 
%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = Nbeam;
mp.P2.compact.Nbeam = mp.P1.compact.Nbeam;
mp.P3.compact.Nbeam = mp.P1.compact.Nbeam;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam;

%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;
mp.NrelayFend = 0; %--How many times to rotate the final image by 180 degrees


%% Optical Layout: Full Model 

%--Pupil Plane Resolutions
mp.P1.full.Nbeam = mp.P1.compact.Nbeam;
mp.P2.full.Nbeam = mp.P1.full.Nbeam;
mp.P3.full.Nbeam = mp.P1.full.Nbeam;
mp.P4.full.Nbeam = mp.P1.full.Nbeam;

%% Pupil Mask Definitions

%--Pupil definition
mp.P1.D = 4; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.

% Both full and compact models
inputs.OD = 1;

% Full model
inputs.Nbeam = mp.P1.full.Nbeam; % number of points across the pupil diameter
inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
mp.P1.full.mask = falco_gen_pupil_Simple(inputs);

% Compact model
inputs.Nbeam = mp.P1.compact.Nbeam; % number of points across usable pupil   
inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam)); % number of points across usable pupil 
mp.P1.compact.mask = falco_gen_pupil_Simple(inputs);

% %--Lyot stop padding
% mp.P4.ODnorm = 183.3529/229.5723; %--Lyot stop OD [Dtelescope]
% inputs.OD = mp.P4.ODnorm; %--Lyot stop OD [Dtelescope]

%% Lyot stop (P4) Definition and Generation

% Inputs common to both the compact and full models
mp.P4.ODnorm = 1;%--Lyot stop OD [Dtelescope]
inputs.OD = mp.P4.ODnorm;%--Lyot stop OD [Dtelescope]

% Full model
inputs.Nbeam = mp.P4.full.Nbeam;
inputs.Npad = 2*round(inputs.Nbeam/2)+2;%2^(nextpow2(3*mp.P4.full.Nbeam));
mp.P4.full.mask = falco_gen_pupil_Simple(inputs); 

% Compact model
inputs.Nbeam = mp.P4.compact.Nbeam;
inputs.Npad = 2*round(inputs.Nbeam/2)+2;%2^(nextpow2(3*mp.P4.compact.Nbeam));
mp.P4.compact.mask = falco_gen_pupil_Simple(inputs); 

clear inputs;

% mp.P4.full.mask = ones(ceil(mp.P4.full.Nbeam+4));
% mp.P4.compact.mask = ones(ceil(mp.P4.compact.Nbeam+4));


%% HLC-Specific Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mp.coro = 'HLC';
% 
% mp.F3.whichMask = 'Seo2019';
% 
% %%-- FPM Material Properties
% 
% mp.F3.substrate = 'fs';
% mp.F3.metal = 'nickel';
% mp.F3.dielectric = 'MgF2'; %'pmgi';
% 
% mp.aoi = 5.5; % Angle of incidence at FPM [deg]
% if(strcmpi(mp.F3.whichMask,'Seo2019'))
%     mp.t_Ti_nm = 3.0; %--Static base layer of titanium beneath any nickel [nm]
% else
%     mp.t_Ti_nm = 0.0; %--Static base layer of titanium beneath any nickel [nm]
% end
% mp.dt_metal_nm = 0.1; %--thickness step size for FPM metal layer (nm)
% mp.t_metal_nm_vec = 0:mp.dt_metal_nm:120; %--nickel thickness range and sampling (nm)
% mp.dt_diel_nm = 0.11;%--thickness step size for FPM dielectric layer  (nm)
% mp.t_diel_nm_vec = 0:mp.dt_diel_nm:120; %--PMGI thickness range and sampling (nm)
% 
% %--Number of waves offset from substrate for reference plane. MUST BE MORE THAN MAX THICKNESS OF THE FPM.
% mp.FPM.d0fac = 4;
% 
% % DM8: FPM Metal Thickness
% if(strcmpi(mp.F3.whichMask,'NiOnly'))
%     mp.dm8.V0coef = 100.5; % Nominal Nickel layer thickness [nm]
% elseif(strcmpi(mp.F3.whichMask,'Seo2019'))
%     mp.dm8.V0coef = 100; % Nominal Nickel layer thickness [nm]
% end
% 
% %--DM8 parameters
% mp.dm8.Vmin = 0;
% mp.dm8.Vmax = 100;
% 
% % DM9: FPM Dielectric thickness
% 
% %--DM9 weights and sensitivities: Used by the controller
% mp.dm9.weight = 4; % Jacobian weight for the FPM dielectric. Smaller weight makes stroke larger by the inverse of this factor.
% mp.dm9.act_sens = 10; %--Change in oomph (E-field sensitivity) of DM9 actuators. Chosen empirically based on how much DM9 actuates during a control step.
% mp.dm9.stepFac = 10;%200; %--Adjust the step size in the Jacobian, then divide back out. Used for helping counteract effect of discretization.
% 
% %--Starting dielectric thicknesses
% if(strcmpi(mp.F3.whichMask,'NiOnly'))
%     mp.dm9.V0coef = 0; % Nominal dielectric layer thickness [nm] 
% elseif(strcmpi(mp.F3.whichMask,'Seo2019'))
%     mp.dm9.V0coef = 0; % Nominal dielectric layer thickness [nm] 
% end
% mp.t_diel_bias_nm = 0; %--Thickness of starting uniform bias layer of PMGI [nm]. % (Requires an outer stop in reality if >0, but will run without it to see if it gives essentially the same result as the EHLC but faster)
% 
% % %--DM9 parameters for 3-fold symmetric Zernikes
% mp.dm9.inf0name = '3foldZern';
% mp.dm9.VtoHavg = 1e-9;
% mp.dm9.maxRadialOrder = 1;

%% FPM (F3) Definition and Generation

mp.F3.full.res = 8;
mp.F3.compact.res = 8; 

mp.coro = 'LC';
mp.F3.whichMask = 'Seo2019';

Fnum = 0.641891/(mp.dm1.Nactbeam*mp.dm1.dm_spacing); % theoretical focal ratio
mp.F3.Rin = 98e-6/(mp.lambda0*Fnum)/2; % radius of inner hard edge of the focal plane mask [lambda0/D]
mp.F3.Rout = Inf;   % radius of outer opaque edge of FPM [lambda0/D]
mp.F3.ang = 180;    % on each side, opening angle [degrees]
mp.FPMampFac = 1.0;%10^(-4.0/2.0); % amplitude transmission of the FPM

mp = dst2_LC_resetFPM(mp);

% Fnum = 1.524/(mp.dm1.Nactbeam*1e-3); % theoretical focal ratio (fOAP=1.524m)
% mp.F3.Rin = 101.6e-6/(mp.lambda0*Fnum)/2;

% mp.F3.Rout = Inf;   % radius of outer opaque edge of FPM [lambda0/D]
% mp.F3.ang = 180;    % on each side, opening angle [degrees]
% mp.FPMampFac = 0.0217;%10^(-3.7/2.0); % amplitude transmission of the FPM
% 
% % Both models
% fpmStruct.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
% fpmStruct.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
% fpmStruct.centering = mp.centering;
% fpmStruct.FPMampFac = mp.FPMampFac; % amplitude transmission of inner FPM spot
% % Full model
% fpmStruct.pixresFPM = mp.F3.full.res;
% mp.F3.full.mask = falco_gen_annular_FPM(fpmStruct);
% % Compact model
% fpmStruct.pixresFPM = mp.F3.compact.res;
% mp.F3.compact.mask = falco_gen_annular_FPM(fpmStruct);

%%

% %--Wavelengths used for Compact Model (and Jacobian Model)
% mp.sbp_weights = ones(mp.Nsbp,1);
% if(strcmpi(mp.estimator,'perfect') && mp.Nsbp>1) %--For design or modeling without estimation: Choose ctrl wvls evenly between endpoints (inclusive) of the total bandpass
%     mp.fracBWsbp = mp.fracBW/(mp.Nsbp-1);
%     mp.sbp_centers = mp.lambda0*linspace( 1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
% %     mp.sbp_weights(1) = 1/2;
% %     mp.sbp_weights(end) = 1/2;
% else %--For cases with estimation: Choose est/ctrl wavelengths to be at subbandpass centers.
%     mp.fracBWsbp = mp.fracBW/mp.Nsbp;
%     mp.fracBWcent2cent = mp.fracBW*(1-1/mp.Nsbp); %--Bandwidth between centers of endpoint subbandpasses.
%     mp.sbp_centers = mp.lambda0*linspace( 1-mp.fracBWcent2cent/2,1+mp.fracBWcent2cent/2,mp.Nsbp); %--Space evenly at the centers of the subbandpasses.
% end

mp.isProbing = false; % initially, not probing 
mp.flagSaveWS = true;
mp.flagSaveEachItr = true;
mp.flagSVD = false; 
mp.full.flagPROPER = false; 
