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


clear all;

close all;

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing with batch process estimation
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT AVAILABLE YET]
% - 'pwp-iekf' for pairwise probing with iterated extended Kalman filter  [NOT AVAILABLE YET]
% mp.estimator = 'pwp-bp';
mp.controller = 'perfect';

%--New variables for estimation:
% - Note: For 360-degree dark hole, must set mp.est.probe.Npairs>=3 and mp.est.probe.axis = 'alternate'.
mp.est.probe.Npairs = 3;%2;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 50;%12;%20;    % Max x/y extent of probed region [actuators].
mp.est.probe.offsetX = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.offsetY = 0;   % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';    % which axis to have the phase discontinuity along [x or y or xy/alt/alternate].
mp.est.probe.gainFudge = 6;%4;%1;     % empirical fudge factor to make average probe amplitude match desired value.

mp.flagSim = true;
mp.layout = 'Fourier';

% mp.dm1.inf_fn = 'influence_BMC_kiloDM_300um_N131.fits';
% mp.dm2.inf_fn = 'influence_BMC_kiloDM_300um_N131.fits';
mp.dm1.inf_fn = 'influence_BMC_kiloDM_300um_N65.fits';
mp.dm2.inf_fn = 'influence_BMC_kiloDM_300um_N65.fits';
% mp.dm1.inf_fn = 'influence_dm5v2.fits';
% mp.dm2.inf_fn = 'influence_dm5v2.fits';

mp.dm1.dm_spacing = 1e-3; %--User defined actuator pitch
mp.dm2.dm_spacing = 1e-3; %--User defined actuator pitch

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';


%% Define Necessary Paths on Your System

%--Library locations
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


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

%%--Coronagraph and Pupil Type
mp.coro = 'Vortex';    %--Tested Options: 'LC','HLC','SPLC','Vortex'
%--Charge of vortex coronagraph as a function of wavelength (can be any
%array of mp.F3.VortexCharge_lambdas and mp.F3.VortexCharge. Will be
%interpolated to find charge for the WFC wavelengths. 
mp.F3.VortexCharge_lambdas = 450e-9:1e-9:550e-9;
mp.F3.VortexCharge = 6*mp.lambda0./mp.F3.VortexCharge_lambdas; 
mp.whichPupil = 'LUVOIR_B_offaxis';%'Simple';
mp.flagApod = false; % Can be simply a sub-aperture 

%--Zernikes to suppress with controller
mp.jac.zerns = 1;%[1 2 3]; %1; %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include at least 1 for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = 2:3; %--Noll indices of Zernikes to compute values for
mp.eval.Rsens = [3, 4;... %--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
                 4, 8];    


%%--Pupil Plane and DM Plane Properties
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)
mp.d_dm1_dm2 = 3; % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
mp.lambda0 = 500e-9; % central wavelength of bandpass (meters)
mp.fracBW = 0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 1; % number of sub-bandpasses across correction band 
mp.Nwpsbp = 1;% number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.


%%--Pupil Masks
switch mp.whichPupil
    case 'Simple' % Can be used to create circular and annular apertures with radial spiders 
        
        mp.P1.D = 4; %--meters, diameter of telescope (This is like HabEx A)
        mp.P1.full.Nbeam = 250; 
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 
        mp.P1.compact.Nbeam = 250;
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex. 
        
        mp.P1.IDnorm = 0;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P1.ODnorm = 1;% Outer diameter (fraction of Nbeam) 
        
        mp.P4.IDnorm = 0;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P4.ODnorm = 0.95;% Outer diameter (fraction of Nbeam) 
        
        mp.P1.Nstrut = 0;% Number of struts 
        mp.P1.angStrut = [];%Array of angles of the radial struts (deg)
        mp.P1.wStrut = []; % Width of the struts (fraction of pupil diam.)
        
        mp.P4.Nstrut = 0;% Number of struts 
        mp.P4.angStrut = [];%Array of angles of the radial struts (deg)
        mp.P4.wStrut = []; % Width of the struts (fraction of pupil diam.)
      
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

        mp.P1.wGap = 6e-3/mp.P1.D; % Fractional gap width
        
end

%% DMs
mp.dm_ind = [1 2]; % vector of which DMs to use for control.
mp.P2.D =     30e-3; % beam diameter at pupil closest to the DMs  (meters)
mp.dm1.Nact = 32; % number of actuators across DM1
mp.dm2.Nact = 32; % number of actuators across DM2
% mp.P2.D =     46.3e-3; % beam diameter at pupil closest to the DMs  (meters)
% mp.dm1.Nact = 48; % number of actuators across DM1
% mp.dm2.Nact = 48; % number of actuators across DM2
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
mp.ctrl.dm9regfacVec = 1;%10.^(-2:1:4);%1/30*10.^(-2:1:2); %--Multiplies with mp.dm9.weight
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


 
%% Final Focal Plane (Fend. Properties
% 
% 
% % %--Specs for Correction (Corr) region and the Scoring (Score) region.
% % mp.Fend.corr.Rin  = 2; %--lambda0/D, inner radius of correction region
% % mp.Fend.score.Rin = 2; %--Needs to be <= that of Correction mask
% % mp.Fend.corr.Rout  = floor(mp.dm1.Nact/2*(1-mp.fracBW/2)); %--lambda0/D, outer radius of correction region
% % mp.Fend.score.Rout = mp.Fend.corr.Rout; %--Needs to be <= that of Correction mask
% % mp.Fend.corr.ang  = 180; %--degrees per side
% % mp.Fend.score.ang = 180; %--degrees per side
% % mp.Fend.sides = 'both'; %--options: 'left', 'right','top','bottom'; any other values produce an annular region 
% 
% 
% % %%--Final Focal Plane (Fend. Properties
% % mp.Fend.res = 3; %--Pixels per lambda_c/D
% % mp.Fend.full.res = 6; %--Pixels per lambda_c/D
% % mp.Fend.FOV = 1 + mp.Fend.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D




%% Part 2: Call the function to define the rest of the variables and initialize the workspace
if(exist('mp','var')==false); mp.dummy = 1; end
mp = falco_config_defaults_VC(mp); %--Load defaults for undefined values

%% Part 3: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Part 4: Run the WFSC trial
out = falco_wfsc_loop(mp);





