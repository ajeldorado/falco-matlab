% Created on 2022-02-23 by Niyati Desai.

%This defaults draft might not be the best, see defaults_VC_simple
%Note both have some parameters commented out because they are filled in with
%SVC mains manually defining them

% Initialize some structures if they don't already exist.
% Generate or load all parameters to run SVC simulation.
%% Misc

%--Record Keeping
mp.SeriesNum = 867;
mp.TrialNum = 5309;

%--Special Computational Settings
mp.flagParfor = false;
mp.useGPU = false;
mp.flagPlot = false;
mp.flagJitter = false;

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
% mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 5;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp-square' for pairwise probing with batch process estimation in a
% square region for one star [original functionality of 'pwp-bp' prior to January 2021]
% - 'pwp-bp' for pairwise probing in the specified rectangular regions for
%    one or more stars
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT TESTED YET]
mp.estimator = 'perfect';

%--New variables for pairwise probing estimation:
mp.est.probe.Npairs = 3;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 12;    % Max x/y extent of probed region [lambda/D].
mp.est.probe.xOffset = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.xiOffset = 0;  %ADDED by Niyati
mp.est.probe.etaOffset = 0; %ADDED by Niyati
mp.est.probe.width = 12; %ADDD by Niyati
mp.est.probe.height = 12; %ADDED by Niyati
mp.est.probe.yOffset = 0;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.

%% Wavefront Control: General

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppress with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = [];%2:3; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = [];%[2,3; 3,4; 4,5]; 

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -6:1:1; %-6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
% % mp.ctrl.dm9regfacVec = 1;        %--Additional regularization factor applied to DM9
   
%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1;
mp.dm2.weight = 0;

%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
mp.controller = 'gridsearchEFC';

% % % GRID SEARCH EFC DEFAULTS     
%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 1; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1;%1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
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
mp.dm1.Nact = 34;%32;               % # of actuators across DM array
mp.dm1.VtoH = 1e-9*ones(mp.dm1.Nact);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = (mp.dm1.Nact/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = (mp.dm1.Nact/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--DM2 parameters
mp.dm2.Nact = 34;%32;               % # of actuators across DM array
mp.dm2.VtoH = 1e-9*ones(mp.dm1.Nact);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;                % clocking of DM surface [degrees]
mp.dm2.xc = (mp.dm2.Nact/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = (mp.dm2.Nact/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm2.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.dm1.Dstop = 100e-3;%mp.dm1.Nact*mp.dm1.dm_spacing;  %--Diameter of iris [meters]
mp.flagDM2stop = true;%false;  %--Whether to apply an iris or not
mp.dm2.Dstop = 0.4*34e-3;%mp.dm1.Nact*mp.dm1.dm_spacing;   %--Diameter of iris [meters]

%--DM separations
mp.d_P2_dm1 = 0;        % distance (along +z axis) from P2 pupil to DM1 [meters]
mp.d_dm1_dm2 = 0.20;   % distance between DM1 and DM2 [meters]


%% Optical Layout: All models

%--Key Optical Layout Choices
mp.flagSim = true;      %--Simulation or not
mp.layout = 'Fourier';  %--Which optical layout to use
mp.coro = 'vortex';

%--Final Focal Plane Properties
mp.Fend.res = 2.5;%3; %--Sampling [ pixels per lambda0/D]
mp.Fend.FOV = 11;%15; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin = 3.0;%2.0;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 3.0;%2.0;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 10;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 180;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'both'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = 0.4*32e-3;%mp.dm1.Nact*mp.dm1.dm_spacing;
mp.P3.D = 0.4*32e-3;%mp.P2.D;
mp.P4.D = 0.4*32e-3;%mp.P2.D;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 600; % INCREASE SAMPLING RES HERE
mp.P2.compact.Nbeam = mp.P1.compact.Nbeam;
mp.P3.compact.Nbeam = mp.P1.compact.Nbeam;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam;  % P4 must be the same as P1 for Vortex. 

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
mp.P1.full.Nbeam = 600; %INCREASE SAMPLING RES HERE
mp.P2.full.Nbeam = mp.P1.full.Nbeam;
mp.P3.full.Nbeam = mp.P1.full.Nbeam;
mp.P4.full.Nbeam = mp.P1.full.Nbeam;  % P4 must be the same as P1 for Vortex. 

% mp.F3.full.res = 6;    % sampling of FPM for full model [pixels per lambda0/D]

%% Entrance Pupil (P1) Definition and Generation

mp.whichPupil =  'Simple';% 'LUVOIR_B'; Used only for run label

%%FROM AVC_defaults
mp.P1.IDnorm = 0.00; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.ODnorm = 1.00;% Outer diameter of the telescope [diameter]
mp.P1.stretch = 1.00; % factor that stretches the horizontal axis to create elliptical beam 
mp.P1.Dfac = 1; %--Factor scaling inscribed OD to circumscribed OD for the telescope pupil.
mp.P1.Nstrut = 0;% Number of struts 
mp.P1.angStrut = [];%Array of angles of the radial struts (deg)
mp.P1.wStrut = []; % Width of the struts (fraction of pupil diam.)

mp.P4.ODnorm = 0.83; %95; %--Lyot stop OD [Dtelescope]
mp.P1.D = 4;%15.2e-3/mp.P4.ODnorm;


%%FROM LC_defaults
% Both full and compact models
inputs.OD = 1.00;

% Full model
inputs.Nbeam = mp.P1.full.Nbeam; % number of points across the pupil diameter
inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
mp.P1.full.mask = falco_gen_pupil_Simple(inputs);

% Compact model
inputs.Nbeam = mp.P1.compact.Nbeam; % number of points across usable pupil   
inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam)); % number of points across usable pupil 
mp.P1.compact.mask = falco_gen_pupil_Simple(inputs);

%COMMENTED THIS PART OUT FOR SIMPLE PUPIL
% mp.P1.D = 7.989; %--circumscribed telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
% 
% %--Generate the entrance pupil aperture
% inputs.centering = mp.centering;
% % Full model:
% inputs.Nbeam = mp.P1.full.Nbeam;
% mp.P1.full.mask = pad_crop(falco_gen_pupil_LUVOIR_B(inputs), 2^(nextpow2(inputs.Nbeam)));
% % Compact model
% inputs.Nbeam = mp.P1.compact.Nbeam;
% mp.P1.compact.mask = pad_crop(falco_gen_pupil_LUVOIR_B(inputs), 2^(nextpow2(inputs.Nbeam)));

%% "Apodizer" (P3) Definition and Generation
mp.flagApod = true;    % false;%--Whether to use an apodizer or not in the FALCO models.

% Inputs common to both the compact and full models
inputs.ID = 0;
inputs.OD =0.84;%1;%0.9;%
inputs.Nstrut = 0;
inputs.angStrut = []; %Angles of the struts 
inputs.wStrut = 0; % spider width (fraction of the pupil diameter)

% Full model only
inputs.Nbeam = mp.P1.full.Nbeam;
inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam)); 
mp.P3.full.mask = falco_gen_pupil_Simple(inputs);

% Compact model only 
inputs.Nbeam = mp.P1.compact.Nbeam;
inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam)); 
mp.P3.compact.mask = falco_gen_pupil_Simple(inputs);


%% Lyot stop (P4) Definition and Generation

% Inputs common to both the compact and full models
inputs.ID = 0;
inputs.OD = mp.P4.ODnorm;

% Full model
inputs.Nbeam = mp.P4.full.Nbeam;
inputs.Npad = 2^(nextpow2(mp.P4.full.Nbeam));
mp.P4.full.mask = falco_gen_pupil_Simple(inputs); 

% Compact model
inputs.Nbeam = mp.P4.compact.Nbeam;
inputs.Npad = 2^(nextpow2(mp.P4.compact.Nbeam));
mp.P4.compact.mask = falco_gen_pupil_Simple(inputs); 


%% VC-Specific Values %%%%%

% mp.F3.VortexCharge = -8; %--Charge of the vortex mask
% mp.F3.phaseMaskType = 'vortex';

% mp.F3.NstepStaircase = 6;
% mp.F3.clocking = 45;
 
% mp.F3.compact.res = 4; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m
% mp.F3.full.res = 16; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m



%For chromaticity (Needs to be cleaned up still)
chromatic = true;

if(chromatic)
    mp.sbp_weights = ones(mp.Nsbp,1);
    if(mp.Nwpsbp==1 && mp.flagSim) %--Set ctrl wavelengths evenly between endpoints (inclusive) of the total bandpass.
        if(mp.Nsbp==1)
            mp.sbp_centers = mp.lambda0;
        else
            mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBW/2, 1+mp.fracBW/2, mp.Nsbp);
            mp.sbp_weights(1) = 1/2; %--Give end sub-bands half weighting
            mp.sbp_weights(end) = 1/2; %--Give end sub-bands half weighting
        end
    else %--For cases with multiple sub-bands: Choose wavelengths to be at subbandpass centers since the wavelength samples will span to the full extent of the sub-bands.
        mp.fracBWcent2cent = mp.fracBW*(1-1/mp.Nsbp); %--Bandwidth between centers of endpoint subbandpasses.
        mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBWcent2cent/2, 1+mp.fracBWcent2cent/2, mp.Nsbp); %--Space evenly at the centers of the subbandpasses.
    end


    mp.F3.phaseScaleFacLambdas = mp.sbp_centers

    mp.F3.phaseScaleFac = mp.lambda0./ mp.F3.phaseScaleFacLambdas
end
