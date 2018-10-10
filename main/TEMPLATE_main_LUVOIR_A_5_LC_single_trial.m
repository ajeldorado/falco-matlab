% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to:
%  1) Specify the key parameter values for a Lyot coronagraph.
%  2) Load the rest of the default settings.
%  3) Save out all the input parameters.
%  4) Run a single trial of WFSC using FALCO.
%
% Modified on 2018-04-05 by A.J. Riggs from HLC to LC.
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
  cd ../ %--Go up one step from the "main" directory to the primary FALCO directory
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

%%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 10; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

mp.controller = 'plannedEFC';  % Controller options: 'EFC' for EFC as an empirical grid search over tuning parametrs, 'conEFC' for constrained EFC using CVX.
% mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

%%--Coronagraph and Pupil Type
mp.coro = 'LC';   %--Tested Options: 'LC','HLC','SPLC','Vortex'
% mp.flagApod = false;
mp.whichPupil = 'LUVOIRA5';

mp.P1.gap_width_m = 0.; %--Width of gaps between primary mirror segments [meters]. Comment out this line to use the default (of 6e-3).
% mp.LS_strut_width = 1.5/100; %--Width of Lyot stop struts [pupil diameters]

%%--Pupil Plane and DM Plane Properties
mp.d_dm1_dm2 = 1; % distance between DM1 and DM2 (meters)
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)

%%--Bandwidth and Wavelength Specs
mp.lambda0 = 600e-9; % central wavelength of bandpass (meters)
mp.fracBW = 0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 1; % number of sub-bandpasses across correction band 
% mp.Nwpsbp = 1;% number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.


%%--Pupil Masks
switch mp.whichPupil

    case{'LUVOIRA5'}  % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        mp.P1.D = 15.2; %14.9760; %--meters, circumscribing diameter of telescope (used only for mas-to-lambda/D conversion)
        mp.P1.Dfac = 15.2/13.7; %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.

        %--Number of pixels across the beam at P1 (the input pupil) and P4 (the Lyot plane)
        mp.P1.full.Nbeam = 400;%500;%800; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P1.compact.Nbeam = 400;%250;%310;%450;
        
        %--Lyot Stop parameters
        mp.P4.IDnorm = 0.40; %--Inner diameter of Lyot stop (relative to circumscribed outer diameter of telescope pupil)
        mp.P4.ODnorm = 0.80; %--Outer diameter of Lyot stop (relative to circumscribed outer diameter of telescope pupil)
        
    case 'LUVOIR_B_offaxis' % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.
        
        mp.P1.full.Nbeam = 800; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for LC. 
        
        mp.P1.compact.Nbeam = 450;
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for LC.

        mp.P3.IDnorm = 0;
        mp.P3.ODnorm = 0.84;
        
        mp.P4.IDnorm = 0; %--Inner diameter of Lyot stop (relative to circumscribed outer diameter of telescope pupil)
        mp.P4.ODnorm = 0.82; %--Outer diameter of Lyot stop (relative to circumscribed outer diameter of telescope pupil)
        
end

%% DMs
mp.dm_ind = [1 2]; % vector of which DMs to use for control.

mp.dm1.Nact = 64; % number of actuators across DM1
mp.dm2.Nact = 64; % number of actuators across DM2
mp.P2.D =     (mp.dm1.Nact-2)*1e-3; % beam diameter at pupil closest to the DMs  (meters). Assumes Xinetics inter-actuator pitch of 1.000 mm.

%--Crop the influence functions used in the DM1 & DM2 Jacobians (but not by PROPER)
inf0 = fitsread('influence_dm5v2.fits');
Nmid = ceil(length(inf0)/2);
Nnew = 41;%51;
infCrop = inf0(Nmid-floor(Nnew/2):Nmid+floor(Nnew/2),Nmid-floor(Nnew/2):Nmid+floor(Nnew/2));
% figure(110); imagesc(log10(abs(inf0)));
% figure(111); imagesc(infCrop);
mp.dm1.inf0 = infCrop;                            
mp.dm2.inf0 = infCrop;   

% %--DM Aperture Mask Properties:
% mp.flagDM1stop = false; %--logical flag whether to include the stop at DM1 or not
% mp.flagDM2stop = false; %--logical flag whether to include the stop at DM2 or not


%% Controller Settings
mp.ctrl.dm9regfacVec = 1;%10.^(-2:1:4);%1/30*10.^(-2:1:2); %--Multiplies with mp.dm_weights(9)
switch mp.controller
    case{'gridsearchEFC'} % 'gridsearchEFC' = empirical grid search over both overall scaling coefficient and log10(regularization)
        % Take images for different log10(regularization) values and overall command gains and pick the value pair that gives the best contrast
        
        mp.dm_ind = [1 2]; %--Which DMs to use and when

        mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
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
            % Column 3: which DMs to use (12, 128, 129, or 1289) for control
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


%% Final Focal Plane (F4) Properties

mp.F4.res = 3; %--Pixels per lambda_c/D
% mp.F4.FOV = 1 + mp.F4.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D



%% Part 2: Call the function to define the rest of the variables and initialize the workspace
if(exist('mp','var')==false); mp.dummy = 1; end
mp = falco_config_defaults_LC(mp); %--Load defaults for undefined values

%% Part 3: Run the WFSC trial
out = falco_wfsc_loop(mp);
falco_wfsc_loop(fn_config,flagPlot); 





