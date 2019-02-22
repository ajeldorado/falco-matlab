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

restoredefaultpath;
rehash toolboxcache;

%% Define Necessary Paths on Your System

pathStem = '~/Repos/';
pathStem2 = '~/Documents/MATLAB/';

% pathStem = '/Users/gruane/Desktop/SCDA/';

%--Library locations
mp.path.falco = [pathStem,'falco-matlab/'];  %--Location of FALCO
mp.path.proper = [pathStem2,'PROPER/'];  %--Location of FALCO; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = [mp.path.falco,'data/brief/']; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = [mp.path.falco,'data/ws/']; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path

%% 

%% Properties for BMC Survey

% %--Variables for survey:
% dz_vec = (3:1.5:10)/10;
% flagWFE_vec = [0, 1]; %--Flag for DM WFE

flagDMwfe = false;
dz = 0.6; % separation distance between DMs [meters]


%--WFE maps and stops on DMs:
mp.dm1.inf_fn = 'influence_BMC_2kDM_400micron_res20.fits';
mp.dm2.inf_fn = 'influence_BMC_2kDM_400micron_res20.fits';

%%--Pupil Masks
mp.P1.full.Nbeam = 937;%309;%250; %--Number of pixel widths across the actual diameter of the beam/aperture (independent of beam centering)
mp.P1.compact.Nbeam = mp.P1.full.Nbeam;

mp.flagDMwfe = flagDMwfe;
if(mp.flagDMwfe)
    mp.dm1.wfe = fitsread(sprintf('~/Data/BMC/wfe/bmc50_dm_wfe_%dpix_pupil.fits',mp.P1.full.Nbeam));
    mp.dm2.wfe = mp.dm1.wfe;
end

mp.dm1.dm_spacing = 400e-6; %--User defined actuator pitch
mp.dm2.dm_spacing = 400e-6; %--User defined actuator pitch

mp.flagDM1stop = false; %--logical flag whether to include the stop at DM1 or not
mp.flagDM2stop = true; %--logical flag whether to include the stop at DM2 or not
mp.dm2.Dstop = 49*mp.dm2.dm_spacing;  %--diameter of circular stop at DM2 and centered on the beam

mp.P2.D =     46.3*mp.dm1.dm_spacing; % beam diameter at pupil closest to the DMs  (meters)
mp.dm1.Nact = 50; % number of actuators across DM1
mp.dm2.Nact = 50; % number of actuators across DM2

%%--Pupil Plane and DM Plane Properties
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)
mp.d_dm1_dm2 = dz; % distance between DM1 and DM2 (meters)


%%--Bandwidth and Wavelength Specs
mp.lambda0 = 550e-9;%575e-9; % central wavelength of bandpass (meters)
mp.fracBW = 0.10;%0.125;%0.10;%0.125;%0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 5; % number of sub-bandpasses across correction band 
mp.Nwpsbp = 1;% number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.
disp(['Wavelengths = ',num2str(mp.lambda0*linspace( 1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp)*1e6)]);



%% Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.useGPU = false; %--whether to use GPUs for Jacobian calculation

mp.flagPlot = true;

%% Step 1: Define any variable values that will overwrite the defaults (in falco_config_defaults_SPLC)

%%--Record Keeping
mp.TrialNum = 1;%8;%6; %--Always use a diffrent Trial # for different calls of FALCO.
mp.SeriesNum = 34; %--Use the same Series # for sets of similar trials.

mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'conEFC' for constrained EFC using CVX. --> DEVELOPMENT ONLY
%mp.controller = 'gridsearchEFC';%--Controller options: 'gridsearchEFC' or 'plannedEFC'
mp.controller = 'plannedEFC';%--Controller options: 'gridsearchEFC' or 'plannedEFC'

% Estimator options:
mp.estimator = 'perfect'; 


%%--Coronagraph and Pupil Type
mp.coro = 'LC';    %--Tested Options: 'LC','HLC','SPLC','Vortex'
mp.flagApod = false;
mp.whichPupil = 'Simple';



%%--Pupil Masks
switch mp.whichPupil
    case 'Simple' % Can be used to create circular and annular apertures with radial spiders 
        
        mp.P1.D = 4; %--meters, diameter of telescope (This is like HabEx A)
%         mp.P1.full.Nbeam = 300; 
%         mp.P1.compact.Nbeam = 300;

        mp.P1.IDnorm = 0;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P1.ODnorm = 1;% Outer diameter (fraction of Nbeam) 
        
        mp.P4.IDnorm = 47.36/227.86;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P4.ODnorm = 156.21/227.86;% Outer diameter (fraction of Nbeam) 
        
        mp.P1.num_strut = 0;% Number of struts 
        mp.P1.strut_angs = [];%Array of angles of the radial struts (deg)
        mp.P1.strut_width = []; % Width of the struts (fraction of pupil diam.)
        
        mp.P4.num_strut = 3;% Number of struts 
        mp.P4.strut_angs = [90 210 330];%Array of angles of the radial struts (deg)
        mp.P4.strut_width = 0.005; % Width of the struts (fraction of pupil diam.)
        
   
end

%% DMs

mp.dm_ind = [1 2 ]; % vector of which DMs to use for control.

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';

%% Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing with batch process estimation
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT AVAILABLE YET]
% - 'pwp-iekf' for pairwise probing with iterated extended Kalman filter  [NOT AVAILABLE YET]
% mp.estimator = 'pwp-bp';
mp.estimator = 'perfect';

%--New variables for estimation:
mp.est.probe.Npairs = 3;%2;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 12;%20;    % Max x/y extent of probed region [actuators].
mp.est.probe.offsetX = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.offsetY = 14;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.

mp.flagSim = true;
mp.layout = 'Fourier';

%% Controller Settings

%--Zernikes to suppress with controller
mp.jac.zerns = 1;%[1 2 3]; %1; %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include at least 1 for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = 2:3; %--Noll indices of Zernikes to compute values for
mp.eval.Rsens = [3, 4;... %--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
                 4, 8];    

mp.ctrl.dm9regfacVec = 1;%10.^(-2:1:4);%1/30*10.^(-2:1:2); %--Multiplies with mp.dm_weights(9)
switch mp.controller
    case{'gridsearchEFC'} % 'EFC' = empirical grid search over both overall scaling coefficient and log10(regularization)
        % Take images for different log10(regularization) values and overall command gains and pick the value pair that gives the best contrast
        
        mp.dm_ind = [1 2]; %--Which DMs to use and when

        mp.ctrl.log10regVec = -5:0.5:-2; %--log10 of the regularization exponents (often called Beta values)
        mp.maxAbsdV = 150;  %--Max +/- delta voltage step for each actuator for DMs 1 and 2
        mp.ctrl.dmfacVec = 1;
        
        %%--WFSC Iterations and Control Matrix Relinearization
        mp.Nitr = 100; %--Number of estimation+control iterations to perform
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
        


%         mp.ctrl.sched_mat = [...
%             repmat([1, 1j, 12, 1, 1], [3, 1, 1]);...
%             repmat([1, -3.5, 12, 0, 0], [2, 1, 1]);...
%             
%             repmat([1, 1j, 12, 1, 1], [3, 1, 1]);...
%             repmat([1, -4.5, 12, 0, 0], [2, 1, 1]);...
%             
%             repmat([1, 1j, 12, 1, 1], [3, 1, 1]);...
%             repmat([1, -5.5, 12, 0, 0], [2, 1, 1]);...
%             
%             repmat([1, 1j, 12, 1, 1], [3, 1, 1]);...
%             repmat([1, -6, 12, 0, 0], [2, 1, 1]);...
%             
%             repmat([1, 1j, 12, 1, 1], [3, 1, 1]);...
%             repmat([1, -6.5, 12, 0, 0], [2, 1, 1]);...
%             
%             repmat([1, 1j, 12, 1, 1], [15, 1, 1]);...
% 
%             ];
        
        mp.ctrl.sched_mat = [...
            repmat([1,1j,12,1,1],[4,1]);...
            repmat([1,1j-1,12,1,1],[25,1]);...
            repmat([1,1j,12,1,1],[2,1]);...
            ];

        [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);
        
        
end
%--Voltage range restrictions
mp.dm1.maxAbsV = 150;
mp.dm2.maxAbsV = 150;

%% Focal plane mask 

mp.F3.Rin = 2.8; % inner hard-edge radius of the focal plane mask, in lambda0/D
mp.FPMampFac = 10^(-3.7);

%% Final Focal Plane (F4) Properties


%--Specs for Correction (Corr) region and the Scoring (Score) region.
mp.F4.corr.Rin  = mp.F3.Rin; %--lambda0/D, inner radius of correction region
mp.F4.score.Rin = mp.F4.corr.Rin; %--Needs to be <= that of Correction mask
mp.F4.corr.Rout  = 10;%floor(mp.dm1.Nact/2*(1-mp.fracBW/2)); %--lambda0/D, outer radius of correction region
mp.F4.score.Rout = mp.F4.corr.Rout; %--Needs to be <= that of Correction mask
mp.F4.corr.ang  = 180; %--degrees per side
mp.F4.score.ang = 180; %--degrees per side
mp.F4.sides = 'both'; %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


%%--Final Focal Plane (F4) Properties
mp.F4.res = 4; %--Pixels per lambda_c/D
mp.F4.FOV = 1 + mp.F4.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D



%% Spatial Weighting of the Control Jacobian
% %--Spatial pixel weighting
% mp.WspatialDef = [mp.F4.corr.Rin, mp.F4.corr.Rin+2, 1];  %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]



%% Part 2: Call the function to define the rest of the variables and initialize the workspace
%if(exist('mp','var')==false); mp.dummy = 1; end
mp = falco_config_defaults_LC(mp); %--Load defaults for undefined values

tag = 'partialTransFPM3_';
mp.runLabel = [tag,mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.F4.corr.Rin),'_OWA',num2str(mp.F4.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Part 3: Run the WFSC trial
out = falco_wfsc_loop(mp);
% falco_wfsc_loop(fn_config,flagPlot); 




