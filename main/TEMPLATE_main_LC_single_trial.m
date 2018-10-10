% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to:
%  1) Specify the key parameter values for a hybrid Lyot coronagraph.
%  2) Load the rest of the default settings.
%  3) Save out all the input parameters.
%  4) Run a single trial of WFSC using FALCO.
%
% Modified on 2018-03-27 by A.J. Riggs from AVC to HLC.
% Modified on 2018-03-22 by A.J. Riggs to have default values that can be
%   overwritten if the variable is already defined in falco_config_defaults_AVC.m.
% Modified on 2018-03-01 by Garreth Ruane and A.J. Riggs to be for a vortex coronagraph.
% Created by A.J. Riggs on 2018-01-08.

%% Go the correct starting directory and add all of FALCO to the Matlab path
clear;

if(~isdeployed)
  pwd0 = fileparts(which(mfilename)); %--Path to this file
  cd(pwd0);
  cd ../
  addpath(genpath(pwd)) %--To find apodizer masks and saved pupils
end

%% Output Data Directories (Change here if desired)
% mp.folders.brief =  %--Config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
% mp.folders.ws = % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%% [OPTIONAL] Start from a previous FALCO trial's DM settings

% fn_prev = 
% temp = load(fn_prev,'out');
% whichItr = 10;
% mp.dm1.V = temp.out.dm1.Vall(:,:,whichItr);
% mp.dm2.V = temp.out.dm2.Vall(:,:,whichItr);
% clear temp

%% Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.useGPU = false; %--whether to use GPUs for Jacobian calculation

mp.flagPlot = true;

%% Step 1: Define any variable values that will overwrite the defaults (in falco_config_defaults_SPLC)

%%--Record Keeping
mp.TrialNum = 1; %--Always use a diffrent Trial # for different calls of FALCO.
mp.SeriesNum = 1; %--Use the same Series # for sets of similar trials.

mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'conEFC' for constrained EFC using CVX. --> DEVELOPMENT ONLY
mp.controller = 'gridsearchEFC'; %'plannedEFC';%--Controller options: 'gridsearchEFC' or 'plannedEFC'

% Estimator options:
mp.estimator = 'perfect'; 


%%--Coronagraph and Pupil Type
mp.coro = 'LC';    %--Tested Options: 'LC','HLC','SPLC','Vortex'
mp.flagApod = false;
mp.whichPupil = 'WFIRST180718';

%%--Pupil Plane and DM Plane Properties
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)
mp.d_dm1_dm2 = 1; % distance between DM1 and DM2 (meters)


%%--Cases for different bandwidth specifications:
% Case 1: Design,Modeling,Testbed. 1 wvl per bandpass. Same wvls in compact and full models
  %--> For testbed, make wvl centers the centers of the sbp instead of
  %linearly from end to end. 
% Case 2: Modeling/Testbed w/ estimation. 1 wvl/sbp in compact model. Nwpsbp wvl/sbp in full model. Need to treat end wavelengths specially in full model.

%%--Bandwidth and Wavelength Specs
% NOTE: Actual wavelengths can be set directly in "mp.sbp_centers" vector in 
% falco_init_ws.m instead if desired. Otherwise, they are auto-generated 
% from the bandwidth and # of wavelengths specified below.
mp.lambda0 = 575e-9; % central wavelength of bandpass (meters)
mp.fracBW = 0.10;%0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 2;  % number of wavelengths or sub-bandpasses (sbp) across entire spectral band 

% %%--Tip/Tilt Control
mp.Ntt = 1; %--Number of tip/tilt locations of the star, including the origin (so always set >=1). 1, 4, or 5
mp.NlamForTT = 0; %--Number of wavelengths to control  tip/tilt at. 0,1, 2, 3, or inf (for all)
mp.TToffset = 1; %--tip/tilt offset (mas)

%%--Pupil Masks
switch mp.whichPupil
    case 'Simple' % Can be used to create circular and annular apertures with radial spiders 
        
        mp.P1.D = 4; %--meters, diameter of telescope (This is like HabEx A)
        mp.P1.full.Nbeam = 500; 
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 
        mp.P1.compact.Nbeam = 350;
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex. 
        
        mp.P1.IDnorm = 0.32;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P1.ODnorm = 1;% Outer diameter (fraction of Nbeam) 
        
        mp.P4.IDnorm = 0.50;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P4.ODnorm = 0.80;% Outer diameter (fraction of Nbeam) 
        
        mp.P1.num_strut = 0;% Number of struts 
        mp.P1.strut_angs = [];%Array of angles of the radial struts (deg)
        mp.P1.strut_width = []; % Width of the struts (fraction of pupil diam.)
        
        mp.P4.num_strut = 0;% Number of struts 
        mp.P4.strut_angs = [];%Array of angles of the radial struts (deg)
        mp.P4.strut_width = []; % Width of the struts (fraction of pupil diam.)
        
    case{'WFIRST180718'}
        mp.P1.D = 2.3631; %--meters, diameter of telescope (used only for mas to lambda/D conversion)
        mp.LS_strut_width = 3.8/100.; % nominal pupil's value is 76mm = 3.216%
        
        mp.P1.full.Nbeam = 324; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P1.compact.Nbeam = 324;
        
        %--Lyot Stop parameters
        mp.P4.IDnorm = 0.50;
        mp.P4.ODnorm = 0.80;
        
    case 'WFIRST_onaxis'
        mp.P1.D = 2.3631; %--meters, diameter of telescope (used only for mas to lambda/D conversion)
        mp.pup_strut_width = 2.61/100;%3.2126/100.; %--0.0261  is nominal from 2014 on-axis (in pupil diameters)
        mp.LS_strut_width = 3.0/100;%3.8/100.; % nominal pupil's value is 76mm = 3.216%

        mp.P1.full.Nbeam = 324; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P1.compact.Nbeam = 324;
        
        %--Lyot Stop parameters
        mp.P4.IDnorm = 0.50;
        mp.P4.ODnorm = 0.80;
      
    case{'LUVOIRA5'}  % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        mp.P1.D = 15.2; %14.9760; %--meters, circumscribing diameter of telescope (used only for mas-to-lambda/D conversion)
        mp.P1.Dfac = 15.2/13.7; %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
        
        mp.P1.full.Nbeam = 500;%1000;%1000; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P1.compact.Nbeam = 250;%400;%300;%400;
        
        %--Lyot Stop parameters
        mp.P4.IDnorm = 0.40;
        mp.P4.ODnorm = 0.80;
        
    case 'LUVOIR_B_offaxis' % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.
        
        mp.P1.full.Nbeam = 250; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 
        
        mp.P1.compact.Nbeam = 250;
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex.

        mp.P3.IDnorm = 0;
        mp.P3.ODnorm = 0.84;
        
        mp.P4.IDnorm = 0;
        mp.P4.ODnorm = 0.82;
        
end

%% DMs
mp.dm_ind = [1 2 ]; % vector of which DMs to use for control.
mp.P2.D =     30e-3; %46.3e-3; % beam diameter at pupil closest to the DMs  (meters)
mp.dm1.Nact = 32;%48; % number of actuators across DM1
mp.dm2.Nact = 32;%48; % number of actuators across DM2

%--Crop the influence functions used in the DM1 & DM2 Jacobians (but not by PROPER)
inf0 = fitsread('influence_dm5v2.fits');
Nmid = ceil(length(inf0)/2);
Nnew = 51;
infCrop = inf0(Nmid-floor(Nnew/2):Nmid+floor(Nnew/2),Nmid-floor(Nnew/2):Nmid+floor(Nnew/2));
% figure(110); imagesc(log10(abs(inf0)));
% figure(111); imagesc(infCrop);
mp.dm1.inf0 = infCrop;                            
mp.dm2.inf0 = infCrop;  

%% Controller Settings
mp.ctrl.dm9regfacVec = 1;%10.^(-2:1:4);%1/30*10.^(-2:1:2); %--Multiplies with mp.dm_weights(9)
switch mp.controller
    case{'gridsearchEFC'} % 'EFC' = empirical grid search over both overall scaling coefficient and log10(regularization)
        % Take images for different log10(regularization) values and overall command gains and pick the value pair that gives the best contrast
        
        mp.dm_ind = [1 2]; %--Which DMs to use and when

        mp.ctrl.log10regVec = -5:1:-2; %--log10 of the regularization exponents (often called Beta values)
        mp.maxAbsdV = 150;  %--Max +/- delta voltage step for each actuator for DMs 1 and 2
        mp.ctrl.dmfacVec = 1;
        
        %%--WFSC Iterations and Control Matrix Relinearization
        mp.Nitr = 20; %--Number of estimation+control iterations to perform
        mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

    case{'plannedEFC'}
        mp.dm_ind = [1 2]; % vector of which DMs to compute Jacobians for at some point (not necessarily all at once or all the time). 

        mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
        mp.maxAbsdV = 150;  %--Max +/- delta voltage step for each actuator for DMs 1 and 2
        mp.ctrl.dmfacVec = 1;
        

        %--CONTROL SCHEDULE. Columns of mp.ctrl.sched_mat are: 
            % Column 1: # of iterations, 
            % Column 2: log10(regularization), 
            % Column 3: which DMs to use (12, 128, 12, or 1289) for control
            % Column 4: flag (0 = false, 1 = true), whether to re-linearize
            %   at that iteration.
            % Column 5: flag (0 = false, 1 = true), whether to perform an
            %   EFC parameter grid search to find the set giving the best
            %   contrast .
            % The imaginary part of the log10(regularization) in column 2 is
            %  replaced for that iteration with the optimal log10(regularization)
            % A row starting with [0, 0, 0, 1...] is for relinearizing only at that time
        
        SetA = ... %--DMs 1 & 2 for 30 iterations. Relinearize every iteration.
            repmat([1, 1j, 12, 1, 1], [20, 1, 1]); 
        SetB = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            10, -3, 12, 0, 0;...
            5,  -4, 12, 0, 0;...
            10, -2, 12, 0, 0;...
            ];
        SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            10, -4, 12, 0, 0;...
            5,  -5, 12, 0, 0;...
            10, -2, 12, 0, 0;...
            ];

        mp.ctrl.sched_mat = [...
            SetA;...
            repmat(SetB,[2,1]);...
            repmat(SetC,[8,1]);...
            ]; 
        
        [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);
        
      
        
    case{'conEFC'} %--Constrained EFC
        mp.dm1.dVpvMax = 30;
        mp.dm2.dVpvMax = 30;
        %mp.dm9.dVpvMax = 40;
        mp.ctrl.dmfacVec = 1;
        mp.ctrl.muVec = 10.^(5); %10.^(8:-1:1);
        
end
%--Voltage range restrictions
mp.dm1.maxAbsV = 150;
mp.dm2.maxAbsV = 150;


%% Final Focal Plane (F4) Properties


% %--Specs for Correction (Corr) region and the Scoring (Score) region.
% mp.F4.corr.Rin  = 2; %--lambda0/D, inner radius of correction region
% mp.F4.score.Rin = 2; %--Needs to be <= that of Correction mask
% mp.F4.corr.Rout  = floor(mp.dm1.Nact/2*(1-mp.fracBW/2)); %--lambda0/D, outer radius of correction region
% mp.F4.score.Rout = mp.F4.corr.Rout; %--Needs to be <= that of Correction mask
% mp.F4.corr.ang  = 180; %--degrees per side
% mp.F4.score.ang = 180; %--degrees per side
% mp.F4.sides = 'both'; %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


% %%--Final Focal Plane (F4) Properties
mp.F4.res = 3; %--Pixels per lambda_c/D
% mp.F4.FOV = 1 + mp.F4.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D


%% Spatial Weighting of the Control Jacobian
% %--Spatial pixel weighting
% mp.WspatialDef = [mp.F4.corr.Rin, mp.F4.corr.Rin+2, 1];  %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]



%% Part 2: Call the function to define the rest of the variables and initialize the workspace
if(exist('mp','var')==false); mp.dummy = 1; end
mp = falco_config_defaults_LC(mp); %--Load defaults for undefined values

%% Part 3: Run the WFSC trial
out = falco_wfsc_loop(mp);
falco_wfsc_loop(fn_config,flagPlot); 




