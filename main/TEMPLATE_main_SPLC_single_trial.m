% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to:
%  1) Specify the key parameter values for a shaped pupil Lyot coronagraph.
%  2) Load the rest of the default settings.
%  3) Save out all the input parameters.
%  4) Run a single trial of WFSC using FALCO.
%
% Modified on 2018-03-27 by A.J. Riggs from AVC to SPLC.
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

%% Step 1: Define any variable values that will overwrite the defaults (in falco_config_defaults_SPLC)

%%--Record Keeping
mp.TrialNum = 1; %--Always use a diffrent Trial # for different calls of FALCO.
mp.SeriesNum = 1; %--Use the same Series # for sets of similar trials.

%%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 10; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.useGPU = false; %--whether to use GPUs for Jacobian calculation

mp.controller = 'EFC';%'conEFC';  % Controller options: 'EFC' for EFC as an empirical grid search over tuning parametrs, 'conEFC' for constrained EFC using CVX.
mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

%%--Coronagraph and Pupil Type
mp.coro = 'SPLC';   %--Tested Options: 'Vortex','LC','SPLC'
mp.whichPupil = 'LUVOIRA5'; %--Tested options: 'WFIRST_onaxis', 'WFIRST20180103','LUVOIRA5'
mp.flagApod = true;
mp.SPname = 'luvoirA5bw10';  %--Shaped pupil identifier. Current options: '32WA194','31WA220', 'luvoirA5bw10'

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

        mp.P3.IDnorm = 0;
        mp.P3.ODnorm = 0.84;
        
        mp.P4.IDnorm = 0;
        mp.P4.ODnorm = 0.82;
        
end

%%--DMs
DM.dm_ind = [1 2]; % vector of which DMs to use for control.
mp.P2.D =     46.3e-3; % beam diameter at pupil closest to the DMs  (meters)
DM.dm1.Nact = 48; % number of actuators across DM1
DM.dm2.Nact = 48; % number of actuators across DM2
DM.dm_weights = ones(9,1);   % vector of relative weighting of DMs for EFC




% %%--Controller Settings
% 
% switch mp.controller
%     case{'EFC'} % 'EFC' = empirical grid search over both overall scaling coefficient and Lagrange multiplier
%         % Take images for different Lagrange multiplier values and overall command gains and pick the value pair that gives the best contrast
%         cp.muVec = 10.^(6:-1:1);
%         cp.dmfacVec = 1;%[0.7, 1]; %[0.5, 1, 2];
%         DM.maxAbsdV = 30; %--Max +/- delta voltage step for each actuator for DMs 1,2, and/or 3
%         
%     case{'conEFC'} %--Constrained EFC, written by He Sun of Princeton University
%         DM.dm1.dVpvMax = 40;
%         DM.dm2.dVpvMax = 40;
%         cp.dmfacVec = 1;
% end
% 
% %--Voltage range restrictions
% DM.dm1.maxAbsV = 250./2.;
% DM.dm2.maxAbsV = 250./2.;
% 
% %%--Tip/Tilt Control
% mp.NlamForTT = 1; %--Number of wavelengths to control  tip/tilt at. 0,1, 2, 3, or inf (for all)
% mp.Ntt = 1; %--Number of tip/tilt offsets, including 0 (so always set >=1). 1, 4, or 5
% mp.TToffset = 1; %--tip/tilt offset (mas)
    





%% Coronagraphic Mask Properties:

% mp.flagDM1stop = false; %--logical flag whether to include the stop at DM1 or not
% mp.flagDM2stop = false; %--logical flag whether to include the stop at DM2 or not


%% Final Focal Plane (F4) Properties


% %--Specs for Correction (Corr) region and the Scoring (Score) region.
% mp.F4.corr.Rin  = 2; %--lambda0/D, inner radius of correction region
% mp.F4.score.Rin = 2; %--Needs to be <= that of Correction mask
% mp.F4.corr.Rout  = floor(DM.dm1.Nact/2*(1-mp.fracBW/2)); %--lambda0/D, outer radius of correction region
% mp.F4.score.Rout = mp.F4.corr.Rout; %--Needs to be <= that of Correction mask
% mp.F4.corr.ang  = 180; %--degrees per side
% mp.F4.score.ang = 180; %--degrees per side
% mp.F4.sides = 'both'; %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


% %%--Final Focal Plane (F4) Properties
% mp.F4.compact.res = 3; %--Pixels per lambda_c/D
% mp.F4.full.res = 6; %--Pixels per lambda_c/D
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

%% Deformable Mirror (DM) Parameters

% %--DM1 parameters
% DM.dm1.Nact = 48; % number of actuators across DM1
% DM.dm1.VtoH = 1*1e-9*ones(DM.dm1.Nact); % Gains: volts to meters in surface height;
% DM.dm1.xtilt = 0;
% DM.dm1.ytilt = 0;
% DM.dm1.zrot = 0; %--clocking angle (degrees)
% DM.dm1.xc = (DM.dm1.Nact/2 - 1/2); % x-center of DM in mm, in actuator widths
% DM.dm1.yc = (DM.dm1.Nact/2 - 1/2); % x-center of DM in mm, in actuator widths
% DM.dm1.edgeBuffer = 1; % Radius (in actuator spacings) outside of pupil to compute influence functions for.
% 
% %--DM2 parameters
% DM.dm2.Nact = 48; % number of actuators across DM1
% DM.dm2.VtoH = 1*1e-9*ones(DM.dm2.Nact); % Gains: volts to meters in surface height;
% DM.dm2.xtilt = 0;
% DM.dm2.ytilt = 0;
% DM.dm2.zrot = 0; %--clocking angle (degrees)
% DM.dm2.xc = DM.dm2.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% DM.dm2.yc = DM.dm2.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% DM.dm2.edgeBuffer = 1; % Radius (in actuator spacings) outside of pupil to compute influence functions for.
% 
% %--DM Actuator characteristics
% DM.dm1.dx_inf0 = 1e-4; % meters, sampling of the influence function;
% DM.dm1.dm_spacing = 1e-3;%0.9906e-3; % meters, pitch of DM actuators
% DM.dm1.inf0 = -1*fitsread('influence_dm5v2.fits');    %  -1*fitsread('inf64.3.fits');                              
% DM.dm2.dx_inf0 = 1e-4; % meters, sampling of the influence function;
% DM.dm2.dm_spacing = 1e-3;%0.9906e-3; % meters, pitch of DM actuators
% DM.dm2.inf0 = -1*fitsread('influence_dm5v2.fits');    







%% Part 2: Call the function to define the rest of the variables and initialize the workspace
if(exist('mp','var')==false); mp.dummy = 1; end
if(exist('cp','var')==false); cp.dummy = 1; end
if(exist('ep','var')==false); ep.dummy = 1; end
if(exist('DM','var')==false); DM.dummy = 1; end

[mp,cp,ep,DM] = falco_config_defaults_SPLC(mp,cp,ep,DM); %--Load values


%% Part 3: Save the config file    
mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(DM.dm_ind)),'DM',num2str(DM.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.F4.corr.Rin),'_OWA',num2str(mp.F4.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];
fn_config = ['data/configs/',mp.runLabel,'.mat'];
save(fn_config)
fprintf('Saved the config file: \t%s\n',fn_config)


%% Part 4: Run the WFSC trial

falco_wfsc_loop(fn_config,'PLOT'); %--Plot progress
% falco_wfsc_loop(fn_config); %--Don't plot anything





