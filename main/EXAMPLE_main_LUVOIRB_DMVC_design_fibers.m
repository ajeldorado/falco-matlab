% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform a DMVC simple design run.
%  1) Load the default model parameters for a vortex.
%  2) Specify the values to overwrite.
%  3) Run a single trial of WFC using FALCO.
%
% REVISION HISTORY:
% --------------
% Modified on 2019-02-26 by A.J. Riggs to load the defaults first.
% ---------------

clear all;
% close all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = 'C:\CoronagraphSims\falco-matlab\';  %--Location of FALCO
mp.path.proper = 'C:\CoronagraphSims\FALCO\lib\PROPER\'; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = 'C:\CoronagraphSims\falco-matlab\data\brief\'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = 'C:\CoronagraphSims\falco-matlab\data\ws\'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path
% rmpath([mp.path.cvx 'lib/narginchk_:']) %--Legend plotting issue if CVX's narginchk function is used instead of Matlab's default function.

%% Step 2: Load default model parameters

EXAMPLE_defaults_LUVOIRB_VC_design_fibers

%% Step 3: Overwrite default values as desired

%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;
mp.flagFiber = true;  %--whether to go place single-mode fibers in the focal plane
mp.flagLenslet = false;  %--whether to go through a lenslet array before using the fibers

%--Record Keeping
mp.SeriesNum = 10;
mp.TrialNum = 4;

mp.lambda0 = 690e-9;
mp.fracBW = 0.40;
mp.Nsbp = 12;
mp.Nitr = 6;
mp.estimator = 'pwp-bp';
mp.est.flagUseJac = false;

mp.fineAlignment_it = 0;

mp.F3.VortexCharge = 8; %--Charge of the vortex mask

%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'Series0009_Trial0004_vortex_LUVOIR_B_offaxis_1DM32_z1_IWA2_OWA15_15lams550nm_BW50_gridsearchEFC_snippet.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

% %--DEBUGGING
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

%--DM settings

mp.dm_ind = [1];

mp.dm1.Nact = 32;
mp.dm1.VtoH = 1*1e-9*ones(mp.dm1.Nact);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = (mp.dm1.Nact/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = (mp.dm1.Nact/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]
% mp.dm1.HminStep = 300e-12; %Expected HCST LSB floor is 300 pm min step size.

mp.dm2.Nact = 32;               % # of actuators across DM array
mp.dm2.VtoH = 1*1e-9*ones(mp.dm2.Nact);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 0;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;                % clocking of DM surface [degrees]
mp.dm2.xc = (mp.dm2.Nact/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = (mp.dm2.Nact/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm2.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]
% mp.dm2.HminStep = mp.dm1.HminStep;

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.dm1.Dstop = 100e-3;  %--Diameter of iris [meters]
mp.flagDM2stop = true;  %--Whether to apply an iris or not
mp.dm2.Dstop = 0.4*50e-3;   %--Diameter of iris [meters]

%--Special settings for fibers

if(mp.flagFiber)
    if(mp.flagLenslet)
        mp.Fend.FOV = 2;
        mp.Fend.res = 15; %Has to be much higher than normal to avoid checkerboarding/ringing when going to F5.
        
        %--Fiber tip plane properties (i.e., focal plane of lenslet(s)
        mp.F5.res = 4;
        mp.F5.FOV = 10;
        mp.F5.fiberPos = [0 0]; %Position of the fiber center in F5 in lambda/D.  
                                %Should be zero unless testing fiber/lenslet misalignments.
        
        %--Lenslet properties
        mp.Fend.lensletWavRad = 1.6; %Radius of the lenslet(s) in lambda_0/D
        mp.Fend.Nlens = 1; %Number of lenslets in Fend
        mp.Fend.x_lenslet = [6];% -3 -3];%[4 9 14 19]; %Lenslet positions in Fend in lambda_0/D
        mp.Fend.y_lenslet = [0];% 5.196 -5.196];%[0 0 0 0];
        mp.lensletFL = 150e-6; %Lenslet focal length in meters
        
        %--Off-axis, incoherent point source (exoplanet)
        mp.x_planet = -mp.Fend.x_lenslet(1); %Position of the exoplanet in lambda_0/D
        mp.y_planet = -mp.Fend.y_lenslet(1);
        %Note that the above coordinates are flipped from the lenslet positions
        %for some damn reason, so input the NEGATIVE of where you want the
        %planet to be.
    else
        mp.Fend.FOV = 20;
        mp.Fend.res = 5;
        
        %--Fiber locations and number
        mp.Fend.Nfiber = 3;
        mp.Fend.x_fiber = [6.1888 -3.0944 -3.0944];%[5.3405 -2.6702 -2.6702]; %Fiber core center positions in lambda_0/D
        mp.Fend.y_fiber = [0 5.3597 -5.3597];%[0 4.625 -4.625];
        
        %--Off-axis, incoherent point source (exoplanet)
        mp.x_planet = mp.Fend.x_fiber(1); %Position of the exoplanet in lambda_0/D
        mp.y_planet = mp.Fend.y_fiber(1);
    end
    
    %--Fiber properties
    mp.fiber.a = 0.507;%0.875;%0.66; %Radius of the fiber core in lambda_0/D
    mp.fiber.a_phys = 1.75e-6; %Physical radius of the fiber core in meters
    mp.fiber.NA = 0.12; %Numerical aperture of the fiber
    
    mp.c_planet = 1e-10; %contrast of exoplanet
    mp.thput_eval_x = mp.x_planet;
    mp.thput_eval_y = mp.y_planet;

end

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);