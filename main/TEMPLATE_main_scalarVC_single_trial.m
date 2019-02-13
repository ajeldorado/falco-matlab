% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to:
%  1) Specify the key parameter values for a vortex coronagraph.
%  2) Load the rest of the default settings.
%  3) Save out all the input parameters.
%  4) Run a single trial of WFSC using FALCO.


clear;

restoredefaultpath;
rehash toolboxcache;

%% Define Necessary Paths on Your System

pathStem = 'D:\Dropbox\Caltech\NASA_project\';

%--Library locations
mp.path.falco = [pathStem,'falco-matlab/'];  %--Location of FALCO
mp.path.proper = [pathStem,'PROPER/'];  %--Location of FALCO; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = [mp.path.falco,'data/brief/']; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = [mp.path.falco,'data/ws/']; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


%% Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path


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

%% Step 1: Define any variable values that will overwrite the defaults (in falco_config_defaults_AVC)

%%--Record Keeping
mp.TrialNum = 1; %--Always use a diffrent Trial # for different calls of FALCO.
mp.SeriesNum = 1; %--Use the same Series # for sets of similar trials.

mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'conEFC' for constrained EFC using CVX. --> DEVELOPMENT ONLY
mp.controller = 'gridsearchEFC';%--Controller options: 'gridsearchEFC' or 'plannedEFC'

%%--Pupil Type
mp.whichPupil = 'Simple';%'Simple';

mp.flagApod = false; % Can be simply a sub-aperture 


%%--Pupil Plane and DM Plane Properties
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)
mp.d_dm1_dm2 = 3; % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
mp.lambda0 = 500e-9; % central wavelength of bandpass (meters)
mp.fracBW = 0.1;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 3; % number of sub-bandpasses across correction band 
mp.Nwpsbp = 1;% number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.

%%--Coronagraph type 
mp.coro = 'Vortex';    %--Tested Options: 'LC','HLC','SPLC','Vortex'

%--Charge of vortex coronagraph. For achromatic vortex, give one value. For
% scalar vortex give the effective charge for each wavelength. 

mp.F3.VortexCharge_lambdas = 475e-9:1e-9:525e-9;
mp.F3.VortexCharge = 6*mp.lambda0./mp.F3.VortexCharge_lambdas; 

%%--Pupil Masks
switch mp.whichPupil
    case 'Simple' % Can be used to create circular and annular apertures with radial spiders 
        
        mp.P1.D = 4; %--meters, diameter of telescope (This is like HabEx A)
        mp.P1.full.Nbeam = 500; 
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 
        mp.P1.compact.Nbeam = 500;
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex. 
        
        mp.P1.IDnorm = 0;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P1.ODnorm = 1;% Outer diameter (fraction of Nbeam) 
        
        mp.P4.IDnorm = 0;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P4.ODnorm = 0.8;% Outer diameter (fraction of Nbeam) 
        
        mp.P1.num_strut = 0;% Number of struts 
        mp.P1.strut_angs = [];%Array of angles of the radial struts (deg)
        mp.P1.strut_width = []; % Width of the struts (fraction of pupil diam.)
        
        mp.P4.num_strut = 0;% Number of struts 
        mp.P4.strut_angs = [];%Array of angles of the radial struts (deg)
        mp.P4.strut_width = []; % Width of the struts (fraction of pupil diam.)
      
    case{'LUVOIRA5'}  % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        mp.P1.D = 15.2; %14.9760; %--meters, circumscribing diameter of telescope (used only for mas-to-lambda/D conversion)
        mp.P1.Dfac = 15.2/13.7; %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
        
        mp.P1.full.Nbeam = 500;%1000;%1000; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P1.compact.Nbeam = 250;%400;%300;%400;
        
    case 'LUVOIR_B_offaxis' % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.
        
        mp.P1.full.Nbeam = 250; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 
        
        mp.P1.compact.Nbeam = 250;
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex.

        mp.P3.apodType = 'Simple';
        
        mp.P3.IDnorm = 0;
        mp.P3.ODnorm = 0.84;
        mp.flagApod = true; % Can be simply a sub-aperture 

        mp.P4.IDnorm = 0;
        mp.P4.ODnorm = 0.82;

        mp.P1.gapWidth = 6e-3/mp.P1.D; % Fractional gap width
        
end

%% DMs
mp.dm_ind = [1 2]; % vector of which DMs to use for control.
mp.P2.D =     46.3e-3; % beam diameter at pupil closest to the DMs  (meters)
mp.dm1.Nact = 48; % number of actuators across DM1
mp.dm2.Nact = 48; % number of actuators across DM2
mp.dm_weights = ones(9,1);   % vector of relative weighting of DMs for EFC

%--Crop the influence functions used in the DM1 & DM2 Jacobians (but not by PROPER)
inf0 = fitsread('influence_dm5v2.fits');
Nmid = ceil(length(inf0)/2);
Nnew = 41;%51;
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

        mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
        mp.maxAbsdV = 150;  %--Max +/- delta voltage step for each actuator for DMs 1 and 2
        mp.ctrl.dmfacVec = 1;
        
        %%--WFSC Iterations and Control Matrix Relinearization
        mp.Nitr = 5; %--Number of estimation+control iterations to perform
        mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

        for Itr = 1:mp.Nitr
            mp.dm_ind_sched{Itr} = mp.dm_ind;
        end
        
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
mp.dm1.maxAbsV = 250;%250./2.;
mp.dm2.maxAbsV = 250;%250./2.;


 
%% Final Focal Plane (F4) Properties
% 
% 
% % %--Specs for Correction (Corr) region and the Scoring (Score) region.
mp.F4.corr.Rin  = 2.5; %--lambda0/D, inner radius of correction region
mp.F4.score.Rin = 3; %--Needs to be <= that of Correction mask
mp.F4.corr.Rout  = 10; %--lambda0/D, outer radius of correction region
mp.F4.score.Rout = mp.F4.corr.Rout; %--Needs to be <= that of Correction mask
% % mp.F4.corr.ang  = 180; %--degrees per side
% % mp.F4.score.ang = 180; %--degrees per side
% % mp.F4.sides = 'both'; %--options: 'left', 'right','top','bottom'; any other values produce an annular region 
% 
% 
% % %%--Final Focal Plane (F4) Properties
% % mp.F4.res = 3; %--Pixels per lambda_c/D
% % mp.F4.full.res = 6; %--Pixels per lambda_c/D
% % mp.F4.FOV = 1 + mp.F4.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D




%% Part 2: Call the function to define the rest of the variables and initialize the workspace
if(exist('mp','var')==false); mp.dummy = 1; end
mp = falco_config_defaults_VC(mp); %--Load defaults for undefined values

%% Part 3: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.F4.corr.Rin),'_OWA',num2str(mp.F4.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Part 4: Run the WFSC trial
out = falco_wfsc_loop(mp);





