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

%% 
LS_ID_vec = 0.45; %(48:1:53)/100;
LS_OD_vec = 0.78; %(78:1:83)/100;

vals_list = allcomb(LS_ID_vec,LS_OD_vec).'; %--dimensions: [2 x length(mp.Nttlam)*length(mp.dm_ind) ]
%     Nvals = size(vals_list,2);

Nsurvey = size(vals_list,2);

for isurvey = 1:Nsurvey

clear mp cp ep DM



mp.P4.IDnorm = vals_list(1,isurvey); %--Lyot stop ID
mp.P4.ODnorm = vals_list(2,isurvey); %--Lyot stop OD

fprintf('******** LS ID = %.2f\t\tLS OD = %.2f\t\t ********\n',mp.P4.IDnorm, mp.P4.ODnorm);


%% Special Computational Settings
mp.flagParfor = true; %true; %--whether to use parfor for Jacobian calculation
mp.useGPU = false; %--whether to use GPUs for Jacobian calculation

mp.flagPlot = true;

%% [OPTIONAL] Start from a previous FALCO trial's DM settings

% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

%% Step 1: Define any variable values that will overwrite the defaults (in falco_config_defaults_SPLC)

%%--Record Keeping
mp.TrialNum = 1;%-1*(0 + isurvey); %--Always use a diffrent Trial # for different calls of FALCO.
mp.SeriesNum = 1; %--Use the same Series # for sets of similar trials.


mp.aoi = 10; % Angle of incidence at FPM [deg]
mp.t_Ti_nm = 3.0; %--Static base layer of titanium beneath any nickel [nm]

mp.dt_metal_nm = 0.25;%1/10; %--thickness step size for FPM metal layer (nm)
mp.t_metal_nm_vec = 0:mp.dt_metal_nm:120; %150; %--nickel thickness range and sampling (nm)
mp.dt_diel_nm = 2/10; %--thickness step size for FPM dielectric layer  (nm)
mp.t_diel_nm_vec = 0:mp.dt_diel_nm:900; %--PMGI thickness range and sampling (nm)

mp.dm8.V0coef = 100;%100/1.42; % Nominal Nickel layer thickness [nm] (the 1.42 disappears when the neighboring actuators are considered)
mp.dm9.V0coef = 390;%240/1.42; % Nominal PMGI layer thickness [nm] (the 1.42 disappears when the neighboring actuators are considered)


% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'conEFC' for constrained EFC using CVX. --> DEVELOPMENT ONLY
mp.controller = 'plannedEFC';%--Controller options: 'gridsearchEFC' or 'plannedEFC'

mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

%%--Coronagraph and Pupil Type
mp.coro = 'HLC';    %--Tested Options: 'LC','HLC','SPLC','Vortex'
mp.flagApod = false;
mp.whichPupil = 'WFIRST180718';%'WFIRST_onaxis';

%%--Pupil Plane and DM Plane Properties
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)
mp.d_dm1_dm2 = 1; % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
mp.lambda0 = 575e-9; % central wavelength of bandpass (meters)
mp.fracBW = 0.10;%0.125;%0.10;%0.125;%0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 6; % number of sub-bandpasses across correction band 
mp.Nwpsbp = 1;% number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.

%--FPM
mp.F3.Rin = 2.7; % maximum radius of inner part of the focal plane mask, in lambda0/D
mp.F3.RinA = mp.F3.Rin; % inner hard-edge radius of the focal plane mask (lambda0/D). Needs to be <= mp.F3.Rin 

%%--Pupil Masks        
%     case{'WFIRST180718'}
mp.P1.full.Nbeam = 250;%350;%250; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
mp.P1.compact.Nbeam = mp.P1.full.Nbeam;

%--Lyot Stop parameters
mp.P4.wStrut = 3.6/100.; % nominal pupil's value is 76mm = 3.216%
% mp.P4.IDnorm = 0.53;%0.50;
% mp.P4.ODnorm = 0.79;%0.80;


%%--Controller Settings
mp.ctrl.dm9regfacVec = 1;%10.^(-2:1:4);%1/30*10.^(-2:1:2); %--Multiplies with mp.dm_weights(9)
switch mp.controller
    case{'gridsearchEFC'} % 'gridsearchEFC' = empirical grid search over both overall scaling coefficient and log10(regularization)
        % Take images for different log10(regularization) values and overall command gains and pick the value pair that gives the best contrast
        
        mp.dm_ind = [1 2 9]; %--Which DMs to use and when

        mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
        mp.maxAbsdV = 150;  %--Max +/- delta voltage step for each actuator for DMs 1 and 2
        mp.ctrl.dmfacVec = 1;
        
        %%--WFSC Iterations and Control Matrix Relinearization
        mp.Nitr = 20; %--Number of estimation+control iterations to perform
        mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

    case{'plannedEFC'}
        mp.dm_ind = [1 2 9]; % vector of which DMs to compute Jacobians for at some point (not necessarily all at once or all the time). 

        mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
        mp.logGmin = -7;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

        mp.maxAbsdV = 150; %100  %--Max +/- delta voltage step for each actuator for DMs 1 and 2
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
        
        SetA = ... %--DMs 1 & 2 for x iterations. Relinearize every iteration.
            repmat([1, 1j, 12, 1, 1], [20, 1]); 
        SetB = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            10, -3, 129, 0, 0;...
            5,  -4, 129, 0, 0;...
            10, -2, 129, 0, 0;...
            ];
        SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            10, -4, 129, 0, 0;...
            5,  -5, 129, 0, 0;...
            10, -2, 129, 0, 0;...
            ];

        mp.ctrl.sched_mat = [...
            SetA;...
            repmat(SetB,[2,1]);...
            repmat(SetC,[8,1]);...
            ]; 
        
        [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);
        
        
        %--ANOTHER OPTION
%         SetA = ... %--DMs 1 & 2 for 30 iterations. Relinearize every iteration.
%             repmat([1, 1j, 12, 1, 1], [30, 1]); 
%         SetB = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
%             [0, 0, 0, 1, 1;...
%             5,  1j, 129, 0, 0;...
%             5,  -2, 129, 0, 0;...
%             5,  -1, 129, 0, 0;...
%             ];
%         SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
%             [0, 0, 0, 1, 1;...
%             5,  1j, 129, 0, 0;...
%             5,  1j-1, 129, 0, 0;...
%             10, -2, 129, 0, 0;...
%             5,  -1, 129, 0, 0;...
%             ];
%         SetD = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
%             [0, 0, 0, 1, 1;...
%             5,  1j+1, 129, 0, 0;...
%             5,  1j, 129, 0, 0;...
%             10, -2, 129, 0, 0;...
%             5,  -1, 129, 0, 0;...
%             ];
% 
%         mp.ctrl.sched_mat = [...
%             SetA;...
%             repmat(SetB,[3,1]);...
%             repmat(SetC,[3,1]);...
%             repmat(SetD,[2,1]);...
%             ];
% 
%         [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

        
    case{'conEFC'} %--Constrained EFC
        mp.dm1.dVpvMax = 30;
        mp.dm2.dVpvMax = 30;
        %mp.dm9.dVpvMax = 40;
        mp.ctrl.dmfacVec = 1;
        mp.ctrl.muVec = 10.^(5); %10.^(8:-1:1);
        
end
%--Voltage range restrictions
mp.dm1.maxAbsV = 250;%250./2.;
mp.dm2.maxAbsV = 250;%250./2.;


% %%--Tip/Tilt Control
% mp.NlamForTT = 1; %--Number of wavelengths to control  tip/tilt at. 0,1, 2, 3, or inf (for all)
% mp.Ntt = 1; %--Number of tip/tilt offsets, including 0 (so always set >=1). 1, 4, or 5
% mp.TToffset = 1; %--tip/tilt offset (mas)



%% DMs
mp.P2.D =     46.3e-3; % beam diameter at pupil closest to the DMs  (meters)
mp.dm1.Nact = 48; % number of actuators across DM1
mp.dm2.Nact = 48; % number of actuators across DM2

mp.dm1.xtilt = 8.986; % degrees
mp.dm2.xtilt = 8.986; % degrees


% mp.flagDM1stop = false; %--logical flag whether to include the stop at DM1 or not
mp.flagDM2stop = true; %--logical flag whether to include the stop at DM2 or not
mp.dm2.Dstop = mp.dm2.Nact*1e-3;  %--diameter of circular stop at DM2 and centered on the beam

%% FPM Properties:



%--DM9 weights and sensitivities
mp.dm_weights = ones(9,1);   % vector of relative weighting of DMs' Jacobians for EFC
mp.dm_weights(9) = 2/3;%1;%10; % Jacobian weight for the FPM dielectric. Smaller weight makes stroke larger by the inverse of this factor.
mp.dm9.act_sens = 10; %--Change in oomph (E-field sensitivity) of DM9 actuators. Chosen empirically based on how much DM9 actuates during a control step.
mp.dm9.stepFac = 200; %--Adjust the step size in the Jacobian, then divide back out. Used for helping counteract effect of discretization.

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

%%--DM9 parameters for 3x3 influence function
mp.dm9.actres = 10;%8;%  number of "actuators" per lambda0/D in the FPM's focal plane. On a square actuator array.
mp.dm9.FPMbuffer = -0.5;%0;%0.2; %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)
mp.dm9.inf0name = '3x3';   % = 1/4*[1, 2, 1; 2, 4, 2; 1, 2, 1];  
%->NOTE: Cannot specify F3 resolution independently of number of DM9
%actuators because there are exactly 2 pixels between neighboring actuator centers.

%% Specs for Correction (Corr) region and the Scoring (Score) region in the final image
mp.F4.corr.Rin  = mp.F3.RinA; %--lambda0/D, inner radius of correction region
% mp.F4.score.Rin = mp.F4.corr.Rin; %--Needs to be >= that of Correction mask
mp.F4.corr.Rout  = 10; %--lambda0/D, outer radius of correction region
% mp.F4.score.Rout = mp.F4.corr.Rout; %--Needs to be <= that of Correction mask
% mp.F4.corr.ang  = 180; %--degrees per side
% mp.F4.score.ang = 180; %--degrees per side
% mp.F4.sides = 'both'; %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


%% Final Focal Plane (F4) Properties
mp.F4.res = 3;%2.5;%3; %--Pixels per lambda_c/D
mp.F4.full.res = 6; %--Pixels per lambda_c/D
% mp.F4.FOV = 1 + mp.F4.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D


%% Tip/Tilt, Spatial, and Chromatic Weighting of the Control Jacobian  #NEWFORTIPTILT
% mp.Ntt = 1; %--Number of tip/tilt offsets, including 0 (so always set >=1). 1, 4, or 5
% mp.NlamForTT = 1; %--Number of wavelengths to compute tip/tilt at. 0,1, 2, 3, or inf (for all)
% mp.TToffset = 1; %--tip/tilt offset in mas
% 
% %--Spatial pixel weighting
% mp.WspatialDef = [mp.F4.corr.Rin, mp.F4.corr.Rin+2, 1];  %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]
% 
% %--Chromatic weighting


% %% Part 2: Call the function to define the rest of the variables and initialize the workspace
% if(exist('mp','var')==false); mp.dummy = 1; end
% if(exist('cp','var')==false); mp.ctrl.dummy = 1; end
% if(exist('ep','var')==false); ep.dummy = 1; end
% if(exist('DM','var')==false); mp.dummy = 1; end
% 
% [mp,cp,ep,DM] = falco_config_defaults_HLC_thinfilm(mp,cp,ep,DM); %--Load defaults for undefined values
% 

% %% Part 3: Save the config file    
% mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
%     mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
%     '_IWA',num2str(mp.F4.corr.Rin),'_OWA',num2str(mp.F4.corr.Rout),...
%     '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
%     '_',mp.controller];
% fn_config = ['data/configs/',mp.runLabel,'.mat'];
% save(fn_config)
% fprintf('Saved the config file: \t%s\n',fn_config)

%% Part 2
mp = falco_config_defaults_HLC(mp);

%% Part 3: Run the WFSC trial

out = falco_wfsc_loop(mp);

% falco_wfsc_loop(fn_config,'PLOT'); %--Plot progress
% falco_wfsc_loop(fn_config); %--Don't plot anything


% fprintf('*** Total time for this trial =     %d seconds    (%.2f minutes) \n\n  ***',toc,toc/60)


end %--END OF SURVEY LOOP



