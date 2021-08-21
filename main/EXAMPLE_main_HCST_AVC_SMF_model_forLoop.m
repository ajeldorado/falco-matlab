% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform a VC simple design run.
%  1) Load the default model parameters for a vortex.
%  2) Specify the values to overwrite.
%  3) Run a single trial of WFC using FALCO.
%
% REVISION HISTORY:
% --------------
% Modified on 2019-04-22 by A.J. Riggs to re-organize the system ID code.
% Modified on 2019-04-12 by He Sun to add the system ID function.
% Modified on 2019-02-26 by A.J. Riggs to load the defaults first.
% ---------------

clear all;

%% Tell Matlab which Python to use from where. Needed for the E-M algorithm.
% setenv('PATH', '/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin') %--Put the parent directory of the python version you want Matlab to use.
% pyversion('/usr/local/opt/python3/bin/python3.6'); %--Put the location of the python version you want Matlab to use.

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = '/Users/jllopsay/Documents/GitHub/falco-matlab/';  %--Location of FALCO
mp.path.proper = '/Users/jllopsay/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '/Users/jllopsay/Documents/GitHub/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '/Users/jllopsay/Documents/GitHub/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
mp.path.mask = '/Users/jllopsay/Documents/GitHub/falco-matlab/lib/masks/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
mp.path.ws_inprogress = '/Users/jllopsay/Documents/GitHub/falco-matlab/data/ws_inprogress/';

% %--Library locations. FALCO and PROPER are required. CVX is optional.
% mp.path.falco = 'C:\Lab\falco-matlab';%'~/Repos/falco-matlab/';  %--Location of FALCO
% mp.path.proper = 'C:\Lab\falco-matlab\proper';%'~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% % mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX
% 
% %%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = 'C:\Lab\falco-matlab\data\configs';%'~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
% mp.path.ws = 'C:\Lab\falco-matlab\data\ws';%'~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path
% rmpath([mp.path.cvx 'lib/narginchk_:']) %--Legend plotting issue if CVX's narginchk function is used instead of Matlab's default function.


%% Step 2: Load default model parameters

EXAMPLE_defaults_HCST_AVC_SMF

mp.flagSaveWS = true;
%% Step 3: Overwrite default values as desired

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 2;

%--WFSC Iterations and Control Matrix Relinearization
mp.controller = 'gridsearchEFC';
mp.Nitr = 45; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1;%1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.dm_ind = [1]; %--Which DMs to use
mp.ctrl.log10regVec = -5:1:2; %--log10 of the regularization exponents (often called Beta values)
mp.use_lastJacStruc = true;

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

% SMF
mp.flagFiber = true;  %--whether to go place single-mode fibers in the focal plane
mp.flagLenslet = false;  %--whether to go through a lenslet array before using the fibers

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0004_Trial0004_vortex_Simple_1DM34_z0.2_IWA3_OWA10_1lams775nm_BW1_gridsearchEFC_360DH_segmented_Iter2of100.mat';
% temp = load(fn_prev);
% mp.dm1.dV = temp.DM1V;
% mp.dm2.dV = temp.DM2V;
% clear temp

%--DEBUGGING IN MONOCHROMATIC LIGHT
mp.fracBW = 0.2;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 11;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
if mp.Nsbp==1
    mp.fracBW = 0.01;
end
%% Other stuff
mp.fineAlignment_it = [];

%--Special settings for fibers

if(mp.flagFiber)
    mp.Fend.FOV = 12;
    mp.Fend.res = 5;
    mp.Fend.Nlens = 1; %Number of lenslets in Fend
    mp.Fend.Nfibers = mp.Fend.Nlens;
    %--Fiber locations and number
    mp.Fend.Nfiber = 1;
    mp.Fend.x_fiber = [-5];%[6.1888 -3.0944 -3.0944];%[5.3405 -2.6702 -2.6702]; %Fiber core center positions in lambda_0/D
    mp.Fend.y_fiber = [0];%[0 5.3597 -5.3597];%[0 4.625 -4.625];
    mp.Fend.x_fiber0 = mp.Fend.x_fiber;%[6.1888 -3.0944 -3.0944];%[5.3405 -2.6702 -2.6702]; %Fiber core center positions in lambda_0/D
    mp.Fend.y_fiber0 = mp.Fend.y_fiber;%[0 5.3597 -5.3597];%[0 4.625 -4.625];

    %--Off-axis, incoherent point source (exoplanet)
    mp.x_planet = mp.Fend.x_fiber(1); %Position of the exoplanet in lambda_0/D
    mp.y_planet = mp.Fend.y_fiber(1);
 
    
    %--Fiber properties
    mp.fiber.a = 0.507;%0.875;%0.66; %Radius of the fiber core in lambda_0/D
    mp.fiber.a_phys = 1.75e-6; %Physical radius of the fiber core in meters
    mp.fiber.NA = 0.12; %Numerical aperture of the fiber
    
    mp.c_planet = 1e-10; %contrast of exoplanet
    mp.thput_eval_x = mp.x_planet;
    mp.thput_eval_y = mp.y_planet;

end

%% Step 4: Generate the label associated with this trial

mp.flagSaveEachItr = false;

mp.flagJitter = true;
jitt_arr = 0.01:0.01:0.07;
numtry = numel(jitt_arr);

for jitt_amp = jitt_arr
    mp.Fend.jitt_amp = jitt_amp;

    label = 'EFCSMF_forJitter';
    mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
        mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
        '_xfiber',num2str(mp.Fend.x_fiber),'_yfiber',num2str(mp.Fend.y_fiber),...
        '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
        '_',mp.controller,'_',num2str(mp.Fend.jitt_amp),'jittAmp_',label];


    %% Step 5: Perform the Wavefront Sensing and Control
    % [out, data_train] = falco_adaptive_wfsc_loop(mp);
    [mp, out] = falco_flesh_out_workspace(mp);
    mp.dm1.V = zeros(mp.dm1.Nact,mp.dm1.Nact);
    [~, out] = falco_wfsc_loop(mp, out);
end
if(mp.flagPlot)
    figure; semilogy(0:mp.Nitr,out.InormHist,'Linewidth',3); grid on; set(gca,'Fontsize',20);
end


