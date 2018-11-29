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


clear all;







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



%% Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
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
mp.SeriesNum = 26; %--Use the same Series # for sets of similar trials.


% mp.aoi = 10; % Angle of incidence at FPM [deg]
% mp.t_Ti_nm = 3.0; %--Static base layer of titanium beneath any nickel [nm]
% 
% mp.dt_metal_nm = 0.25;%1/10; %--thickness step size for FPM metal layer (nm)
% mp.t_metal_nm_vec = 0:mp.dt_metal_nm:120; %150; %--nickel thickness range and sampling (nm)
% mp.dt_diel_nm = 2/10; %--thickness step size for FPM dielectric layer  (nm)
% mp.t_diel_nm_vec = 0:mp.dt_diel_nm:900; %--PMGI thickness range and sampling (nm)
% 
% mp.dm8.V0coef = 100;%100/1.42; % Nominal Nickel layer thickness [nm] (the 1.42 disappears when the neighboring actuators are considered)
% mp.dm9.V0coef = 390;%240/1.42; % Nominal PMGI layer thickness [nm] (the 1.42 disappears when the neighboring actuators are considered)

%--Note: dm8 amplitude goes from zero to one, but for simplicity with the
%surface it is used as 1-dm8amplitude for the FPM
mp.dm8.Vmin = 0; % 1-mp.dm8.Vmin is the maximum allowed FPM amplitude
mp.dm8.Vmax = 1 - 1e-4; % 1-mp.dm8.Vmax is the minimum allowed FPM amplitude

mp.dm8.V0coef = 1 - 1e-2; % starting occulting spot amplitude amplitude
mp.dm9.V0coef = -(147/360)*575; %-100; % starting occulting spot phase (uniform) [nm of phase]


% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'conEFC' for constrained EFC using CVX. --> DEVELOPMENT ONLY
% mp.controller = 'SM-AMPL';%--Controller options: 'gridsearchEFC' or 'plannedEFC'
mp.controller = 'plannedEFC';%--Controller options: 'gridsearchEFC' or 'plannedEFC'

mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

%%--Coronagraph and Pupil Type
mp.coro = 'FOHLC';    %--Tested Options: 'LC','HLC','SPLC','Vortex','FOHLC'
mp.flagApod = false;
mp.whichPupil = 'WFIRST180718';%

%%--Pupil Plane and DM Plane Properties
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)
mp.d_dm1_dm2 = 1; % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
mp.lambda0 = 575e-9; % central wavelength of bandpass (meters)
mp.fracBW = 0.10;%0.125;%0.10;%0.125;%0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 1;%5; % number of sub-bandpasses across correction band 
mp.Nwpsbp = 1;% number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.

%--FPM
mp.F3.Rin = 2.7;%4;%2.7; % maximum radius of inner part of the focal plane mask, in lambda0/D
mp.F3.RinA = 2.7;%mp.F3.Rin; % inner hard-edge radius of the focal plane mask (lambda0/D). Needs to be <= mp.F3.Rin 

%%--Pupil Masks        
mp.P1.full.Nbeam = 250;%350;%250; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
mp.P1.compact.Nbeam = mp.P1.full.Nbeam;

%--Lyot Stop parameters
mp.P4.wStrut = 3.6/100.; % nominal pupil's value is 76mm = 3.216%
mp.P4.IDnorm = 0.45; %--Lyot stop ID
mp.P4.ODnorm = 0.78; %--Lyot stop OD

%%--Controller Settings
mp.ctrl.dm9regfacVec = 1;%10.^(-2:1:4);%1/30*10.^(-2:1:2); %--Multiplies with mp.dm_weights(9)
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

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
        mp.dm_ind = [1 2 8 9];%[1 2 8];%8;%[1 2 8 9]; % vector of which DMs to compute Jacobians for at some point (not necessarily all at once or all the time). 

        mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)

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
        SetD = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
           [0, 0, 0, 1, 0;...
           10, -5, 129, 0, 0;...
           5,  -6, 129, 0, 0;...
           10, -2, 129, 0, 0;...
           ];

%         mp.ctrl.sched_mat = [...
%             SetA;...
%             repmat(SetB,[2,1]);...
%             repmat(SetC,[8,1]);...
%             ]; 
        
       SetA2 = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
       SetB2 = [1, -5, 12, 1, 0];
       SetC2 = [1, 1j, 12, 1, 1];
       
       mp.ctrl.sched_mat = [...
           repmat(SetA2,[5,1]);...
           repmat(SetB2,[3,1]);...
           repmat(SetC2,[3,1]);...
           repmat(SetB,[2,1]);...
           repmat(SetC,[4,1]);...
           repmat(SetD,[4,1]);...
           ];
       
%         mp.ctrl.sched_mat = [...
%             repmat([1, 1j, 128, 1, 1], [3, 1]);
%             ]; 
        
        
        [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);
        
        
    case{'SM-AMPL'} %--Constrained EFC with AMPL
        
        addpath('~/bin/amplapi/matlab/')
        
        mp.dm_ind = 8;%[1 2 8 9]; %--Which DMs to use and when
        mp.ctrl.log10regVec = -5:1:-2;%  -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
        mp.ctrl.dmfacVec = 1;

        %--Delta voltage range restrictions
        mp.dm1.maxAbsdV = 50;%30;
        mp.dm2.maxAbsdV = 50;%30;
        mp.dm9.maxAbsdV = 50;%40;
        
        %--Absolute voltage range restrictions
        mp.dm1.maxAbsV = 150;%250;
        mp.dm2.maxAbsV = 150;%250;
        
    case{'conEFC'} %--Constrained EFC
        mp.dm1.dVpvMax = 30;
        mp.dm2.dVpvMax = 30;
        mp.dm9.dVpvMax = 40;
        mp.ctrl.dmfacVec = 1;
        mp.ctrl.muVec = 10.^(5); %10.^(8:-1:1);
        
end


%--Delta voltage range restrictions
mp.dm8.maxAbsdV = 0.2;

%--Voltage range restrictions
mp.dm1.maxAbsV = 1000;
mp.dm2.maxAbsV = 1000;

%%

mp.dm_ind = 8;%[1 8];



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
mp.dm_weights(9) = 1;%2/3;%1;%10; % Jacobian weight for the FPM phase. Smaller weight makes stroke larger by the inverse of this factor.
% % mp.dm9.act_sens = 10; %--Change in oomph (E-field sensitivity) of DM9 actuators. Chosen empirically based on how much DM9 actuates during a control step.
% % mp.dm9.stepFac = 200; %--Adjust the step size in the Jacobian, then divide back out. Used for helping counteract effect of discretization.

%--DM8 weight
mp.dm_weights(8) = 1e-2;%0.1;%1e-3 % Jacobian weight for the FPM amplitude. Smaller weight makes stroke larger by the inverse of this factor.
mp.dm8.act_sens = 1;


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
mp.dm9.actres = 5;%10;%8;%  number of "actuators" per lambda0/D in the FPM's focal plane. On a square actuator array.
mp.dm9.FPMbuffer = 0;%-0.5;%0;%0.2; %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)
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
% mp.F4.FOV = 1 + mp.F4.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D


%% Tip/Tilt, Spatial, and Chromatic Weighting of the Control Jacobian  #NEWFORTIPTILT
% 
% %--Spatial pixel weighting
% mp.WspatialDef = [mp.F4.corr.Rin, mp.F4.corr.Rin+2, 1];  %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%% Part 2
mp = falco_config_defaults_FOHLC(mp);

%% Part 3: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.F4.corr.Rin),'_OWA',num2str(mp.F4.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

% %% Part 4: Run the WFSC trial
% out = falco_wfsc_loop(mp);


%% Part 4: Initialization function (copied from falco_wfsc_loop.m)

%--Sort out file paths and save the config file    

%--Add the slash or backslash to the FALCO path if it isn't there.
if( strcmp(mp.path.falco(end),'/')==false || strcmp(mp.path.falco(end),'\')==false )
    mp.path.falco = [mp.path.falco filesep];
end

mp.path.dummy = 1; %--Initialize the folders structure in case it doesn't already exist

%--Store minimal data to re-construct the data from the run: the config files and "out" structure after a trial go here
if(isfield(mp.path,'config')==false)
    mp.path.config = [mp.path.falco filesep 'data' filesep 'config' filesep];     
end

%--Entire final workspace from FALCO gets saved here.
if(isfield(mp.path,'ws')==false)
    mp.path.ws = [mp.path.falco filesep 'data' filesep 'ws' filesep];      
end

%--Save the config file
fn_config = [mp.path.config mp.runLabel,'_config.mat'];
save(fn_config)
fprintf('Saved the config file: \t%s\n',fn_config)


%--Get configuration data from a function file
[mp,out] = falco_init_ws(fn_config);

%%
%--Update the number of elements used per DM
if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); end
if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); end
if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); end
if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); end

%% Part 5: Compute the Jacobian using model_Jacobian

modvar.flagCalcJac = true; 
modvar.wpsbpIndex = mp.wi_ref;
modvar.whichSource = 'star'; 

%--Re-initialize the Jacobian arrays to full size
G1=zeros(1,1,mp.jac.Nmode); G2=zeros(1,1,mp.jac.Nmode); G3=zeros(1,1,mp.jac.Nmode); G4=zeros(1,1,mp.jac.Nmode); G5=zeros(1,1,mp.jac.Nmode); G6=zeros(1,1,mp.jac.Nmode); G7=zeros(1,1,mp.jac.Nmode); G8=zeros(1,1,mp.jac.Nmode); G9=zeros(1,1,mp.jac.Nmode); %--Initialize for bookkeeping in cells later 
if(any(mp.dm_ind==1)); G1 = zeros(mp.F4.corr.Npix,mp.dm1.NactTotal,mp.jac.Nmode); end % control Jacobian for DM1
if(any(mp.dm_ind==2)); G2 = zeros(mp.F4.corr.Npix,mp.dm2.NactTotal,mp.jac.Nmode); end % control Jacobian for DM2
if(any(mp.dm_ind==8)); G8 = zeros(mp.F4.corr.Npix,mp.dm8.NactTotal,mp.jac.Nmode); end % control Jacobian for DM8
if(any(mp.dm_ind==9)); G9 = zeros(mp.F4.corr.Npix,mp.dm9.NactTotal,mp.jac.Nmode); end % control Jacobian for DM9
    
%--Compute the number of total actuators for all DMs used. 

GallCell1 = {squeeze(G1(:,:,1)),squeeze(G2(:,:,1)),squeeze(G3(:,:,1)),squeeze(G4(:,:,1)),squeeze(G5(:,:,1)),squeeze(G6(:,:,1)),squeeze(G7(:,:,1)),squeeze(G8(:,:,1)),squeeze(G9(:,:,1))}; % Create the cell array. Placeholders for non-existent Jacobians to have consistent numbering
NeleAll = 0;
NeleVec = []; %--Vector of total number of used actuators for each used DM
for ii=1:numel(mp.dm_ind)
    dm_index = mp.dm_ind(ii);
    NeleAll = NeleAll + size(GallCell1{dm_index},2);
    NeleVec = [NeleVec; size(GallCell1{dm_index},2) ];
end
clear GallCell1 %--Save RAM

%--Compute the control Jacobians for each DM
jacStruct =  model_Jacobian(mp);
if(any(mp.dm_ind==1)); G1_jac = jacStruct.G1; end
if(any(mp.dm_ind==2)); G2_jac = jacStruct.G2; end
if(any(mp.dm_ind==8)); G8_jac = jacStruct.G8; end
if(any(mp.dm_ind==9)); G9_jac = jacStruct.G9; end
clear jacStruct  %--Save RAM



%% Part 6: Compute the Jacobian using model_compact and differencing

Vfrac = 1e-4; %--Want delta voltage to be tiny for DM1 and DM2 to stay linear

%--Re-initialize the Jacobian arrays to full size
if(any(mp.dm_ind==1)); G1 = zeros(mp.F4.corr.Npix,mp.dm1.NactTotal,mp.jac.Nmode); end % control Jacobian for DM1
if(any(mp.dm_ind==2)); G2 = zeros(mp.F4.corr.Npix,mp.dm2.NactTotal,mp.jac.Nmode); end % control Jacobian for DM2
if(any(mp.dm_ind==8)); G8 = zeros(mp.F4.corr.Npix,mp.dm8.NactTotal,mp.jac.Nmode); end % control Jacobian for DM8
if(any(mp.dm_ind==9)); G9 = zeros(mp.F4.corr.Npix,mp.dm9.NactTotal,mp.jac.Nmode); end % control Jacobian for DM9
    



for tsi=1:mp.jac.Nmode
    
    Ein = ones(mp.P1.compact.Narr);
    
    modvar.wpsbpIndex = 0; %--Dummy index since not needed in compact model
    
    modvar.sbpIndex = 1;% 
%     modvar.ttIndex = mp.Wttlam_ti(tsi);
%     modvar.flagCalcJac = 0; 
    modvar.whichSource = 'star';     
    
    normFac = mp.F4.compact.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
    flagEval = false;             % flag to use a different (usually higher) resolution at final focal plane for evaluation
    
    lambda = mp.sbp_centers(modvar.sbpIndex);
    
    Eunpoked = model_compact(mp, modvar);
    EunpokedVec = Eunpoked(mp.F4.corr.inds);
    
    
    
    %--DM1
    fprintf('Starting Jacobian calculation with compact model for DM1...'); tic
    if(any(mp.dm_ind==1))
        whichDM = 1;
        parfor iact=1:mp.dm1.NactTotal
            %G1(:,iact,tsi) = func_validate_Jacobian_with_compact_model(iact,whichDM,Vfrac,EunpokedVec,mp, DM, modvar);
            G1(:,iact,tsi) = func_validate_Jacobian_with_compact_model_FOHLC(iact,whichDM, Vfrac, EunpokedVec, mp,  lambda, normFac, Ein)
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)
%     
%     %--DM2
%     fprintf('Starting Jacobian calculation with compact model for DM2...'); tic
%     if(any(mp.dm_ind==2))
%         whichDM = 2;
%         parfor iact=1:mp.dm2.NactTotal
%             G2(:,iact,tsi) = func_validate_Jacobian_with_compact_model(iact,whichDM,Vfrac,EunpokedVec,mp, DM, modvar);
%         end
%     end
%     fprintf('done. Time = %.1f sec.\n',toc)


    %--DM8
    fprintf('Starting Jacobian calculation with compact model for DM8...'); tic
    if(any(mp.dm_ind==8))
        whichDM = 8;
        VfracDM9 = 1;
        parfor iact=1:mp.dm8.NactTotal
           G8(:,iact,tsi) = func_validate_Jacobian_with_compact_model_FOHLC(iact,whichDM, Vfrac, EunpokedVec, mp,  lambda, normFac, Ein)
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)


    %--DM9
    fprintf('Starting Jacobian calculation with compact model for DM9...'); tic
    if(any(mp.dm_ind==9))
        whichDM = 9;
        VfracDM9 = 10;
        parfor iact=1:mp.dm9.NactTotal
           % G9(:,iact,tsi) = func_validate_Jacobian_with_compact_model(iact,whichDM,VfracDM9,EunpokedVec,mp, DM, modvar);
           G9(:,iact,tsi) = func_validate_Jacobian_with_compact_model_FOHLC(iact,whichDM, Vfrac, EunpokedVec, mp,  lambda, normFac, Ein)
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)
        
end

if(any(mp.dm_ind==1));  G1_compact = G1; clear G1;  end
if(any(mp.dm_ind==8));  G8_compact = G8; clear G8;  end
if(any(mp.dm_ind==9));  G9_compact = G9; clear G9;  end

%%
%% DM1 Comparison, Compact Model
%--Compare overall Jacobian
figure(1); imagesc(abs(G1_jac)); colorbar; set(gca,'Fontsize',20);
figure(2); imagesc(abs(G1_compact)); colorbar; set(gca,'Fontsize',20);
figure(3); imagesc(abs(G1_jac-G1_compact)); colorbar; set(gca,'Fontsize',20);
figure(4); imagesc(abs(G1_jac-G1_compact)/max(abs(G1_compact(:)))); colorbar; set(gca,'Fontsize',20); %--Normalized error
% figure(5); imagesc(abs(G1_jac-G1_compact)./abs(G1_compact),[0 1]); colorbar; set(gca,'Fontsize',20);

%--Compare Jacobian for just one actuator
whichAct = 1024;

Etemp = model_compact(mp, modvar);
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.F4.corr.inds) = G1_jac(:,whichAct,1);
Etemp2(mp.F4.corr.inds) = G1_compact(:,whichAct,1);


figure(11); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(12); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(13); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(14); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error
 
% figure(15); imagesc(angle(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(real(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(imag(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);


figure(18); imagesc(angle(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(19); imagesc(angle(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);

% %--Reset the voltage maps to the starting point.
% mp.dm1.V = DM1V0;
% mp.dm2.V = DM2V0;
% mp.dm8.V = DM8V0;
% mp.dm9.V = DM9V0;


%% DM8 Comparison, Compact Model

jac2D = zeros(sqrt(mp.dm8.NactTotal));
jac2D(:) = sum(abs(G8_jac).^2,1);
figure(20); imagesc(abs(jac2D)); colorbar; set(gca,'Fontsize',20);

jac2Db = zeros(sqrt(mp.dm8.NactTotal));
jac2Db(:) = sum(abs(G8_compact).^2,1);
figure(30); imagesc(abs(jac2Db)); colorbar; set(gca,'Fontsize',20);


%--Compare Jacobian for just one actuator
whichAct = 14*28+14;%13*54+50;%27+54*27;%1024;

Etemp = model_compact(mp, modvar); %--Just used to get the right sized image
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.F4.corr.inds) = G8_jac(:,whichAct,1);
Etemp2(mp.F4.corr.inds) = G8_compact(:,whichAct,1);


figure(31); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(32); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(33); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(34); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error
 
% figure(15); imagesc(angle(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(real(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(imag(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);


figure(38); imagesc(angle(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(39); imagesc(angle(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);


%% DM9 Comparison, Compact Model
% 
% jac2D = zeros(sqrt(mp.dm9.NactTotal));
% jac2D(:) = sum(abs(G9_jac).^2,1);
% figure(20); imagesc(abs(jac2D)); colorbar; set(gca,'Fontsize',20);
% 
% jac2Db = zeros(sqrt(mp.dm9.NactTotal));
% jac2Db(:) = sum(abs(G9_compact).^2,1);
% figure(30); imagesc(abs(jac2Db)); colorbar; set(gca,'Fontsize',20);


%--Compare Jacobian for just one actuator
whichAct = 13*54+50;%27+54*27;%1024;

Etemp = model_compact(mp, DM, modvar); %--Just used to get the right sized image
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.F4.corr.inds) = G9_jac(:,whichAct,1);
Etemp2(mp.F4.corr.inds) = G9_compact(:,whichAct,1);


figure(31); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(32); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(33); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(34); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error
 
% figure(15); imagesc(angle(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(real(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(imag(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);


figure(38); imagesc(angle(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(39); imagesc(angle(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
