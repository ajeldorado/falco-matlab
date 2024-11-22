% Copyright 2018 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to run WFSC with a vortex using implicit EFC (IEFC) and a Fourier basis set.

% clear

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
addpath(genpath('/home/hcst/falco-matlab'));% savepath;
addpath('/home/hcst/HCST_scripts/EFC');
%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false

clearvars -except bench 

flag_loop = false;
flag_efc360is=false;
flag_smallDH4smf=false;
flag_apodizer_avc=false;

%************
lambda0 = 780*1e-9;
%************
%% Step 2: Load default model parameters

mp.flagFiber = true;
EXAMPLE_defaults_HCST_VVC 
mp.flagSim = false;      %--Simulation or not
mp.testbed = 'HCST';

mp.Nitr = 55;

%% SMF Params
mp.flagLenslet = false;

mp.Fend.x_fiber = [bench.FIUstages.smf_angular_pos_efc1(1)];%[5.3405 -2.6702 -2.6702]; %Fiber core center positions in lambda_0/D
mp.Fend.y_fiber = [bench.FIUstages.smf_angular_pos_efc1(2)];%[0 4.625 -4.625];
mp.Fend.Nfiber = numel(mp.Fend.x_fiber);

mp.fiber.a = 0.507;%0.875;%0.66; %Radius of the fiber core in lambda_0/D
mp.fiber.a_phys = 1.75e-6; %Physical radius of the fiber core in meters
mp.fiber.NA = 0.12; %Numerical aperture of the fiber

%% HCST preparation
% Flags; all flags false --> regular EFC run
flag_svc = false;
efc_imageSharpening = false;
flagFieldStop = true;
mp.flag_lc = false;

% path2flatMap =[bench.info.HCST_DATA_DIR,'is_fiu_results/2023Jun29/'];
% path2flatMap = [bench.info.HCST_DATA_DIR,'zygo_flat/2019Jun26/'];
path2flatMap = [bench.info.HCST_DATA_DIR,'is_results/2023Sep15/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'is_fiu_results/2023Jan12/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'efcsmf_results/2023Sep09/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'efc_results/2023Sep18/']; %!!!!
path2startingMap =[bench.info.HCST_DATA_DIR,'efc_results/2023Oct19/'];
path2startingMap =[bench.info.HCST_DATA_DIR,'efc_results/2023Dec05/'];

NDfilter_FWpos = 2;
FWpos = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsbp = 3; % Number of sub-bandpasses, AKA number of wavelengths in full bandpass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FWpos~=NDfilter_FWpos
    NDfilter_cal = 28.6082 * 4.9367 * 4.4753 * 4.6956;%23.16; %>1
%     NDfilter_cal = 28.6082 * 4.9367 * 4.4753;% * 4.6956;%23.16; %>1
else
    NDfilter_cal = 1;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tint_offAxis = 2.5e-2*ones(Nsbp,1); %2.5e-2*ones(Nsbp,1); %; % 1e-3; %Integration time for off-axis PSF
mp.tint_efc = 10e-0*ones(Nsbp,1); %30; %1e-1; %30; %35e-0;%1e-1;%[2,1,1,0.5,1]*0.05; % for the control and evalution
mp.tint_est = 5e-1*ones(Nsbp,1); %2.5e-0*ones(Nsbp,1);%5e-1; % time of integration for the estimation/sensing
% tint_offAxis(1:2) = 1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bench.info.sbp_texp = mp.tint_efc;
mp.tint = mp.tint_efc;
mp.peakPSFtint = tint_offAxis;
mp.NDfilter_cal = NDfilter_cal;

% Source parameters
bench.info.source =  'nkt';%'laser';%

if Nsbp>1
    mp.fracBW = 0.1;  % bandwidth of experiment
    mp.sbp_centers = lambda0*linspace( 1-mp.fracBW/2,1+mp.fracBW/2,Nsbp);
    bench.info.sbp_width = 14e-9*ones(Nsbp,1);%[14,14,14,14,14]*1e-9;% 3e-9;%-  %--Width of each sub-bandpass on testbed (meters)
%     bench.info.sbp_width([4:6]) = 6e-9; %
        gainfactor = 0.0135; % Arielle!
        mp.dm1.VtoH = (4e-7*ones(mp.dm1.Nact) * 2*sqrt(2)).*gainmap*gainfactor;%
        mp.est.probe.gainFudge = 0.1;
else
    mp.fracBW = 0.01; % bandwidth of experiment; when Nsbp is 1 --> narrowband or monochromatic  
    mp.sbp_centers = lambda0;
    bench.info.sbp_width = [14]*1e-9; % bandpass of the sub-bandpass [nm]
    if mp.flagFiber
%         gainfactor = 0.005;
        gainfactor = 0.01;
        mp.dm1.VtoH = (4e-7*ones(mp.dm1.Nact) * 2*sqrt(2)).*gainmap*gainfactor;%
        mp.est.probe.gainFudge = 0.1;
    else
        mp.est.probe.gainFudge = 1;
    end
end

frameSize = 180;

datelabel = datestr(now,'yyyymmmddTHHMM');
label_progress = [num2str(Nsbp),'lams',num2str(round(1e9*lambda0)),'nm_BW',num2str(mp.fracBW*100),'_VVC_',datelabel];

falcoPrepare_smfv2

mp.bench = bench;
% bench.info.sbp_width = sbp_width ; %--Width of each sub-bandpass on testbed (meters)
bench.info.sbp_texp = mp.tint_efc;% Exposure time for each sub-bandpass (seconds)
bench.info.PSFpeaks = PSFpeaks;% counts per second 

%% Step 3: Overwrite default values as desired

%--Whether to perform a model-based (instead of empirical) grid search for the controller
mp.ctrl.flagUseModel = false; 

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;


%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;
mp.runLabel = sprintf('Series%04d_Trial%04d', mp.SeriesNum, mp.TrialNum);

% %--Use just 1 wavelength for initial debugging/testing of code
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
% mp.Nitr = 3; %--Number of wavefront control iterations

% mp.dm1.fourier_spacing = 1.0; %0.5; % Center-to-center spacing between Fourier modes in the focal plane. [lambda/D]
% mp.dm1.fourier_gridType = 'hex';  % Options: 'hex' or 'square'. 'hex' has a denser packing
% xiMin = mp.Fend.corr.Rin-1;
% clocking = mp.Fend.clockAngDeg - 90; % -90 for 'bottom' dark hole.
% [mp.dm1.fourier_basis_xis , mp.dm1.fourier_basis_etas] = falco_choose_fourier_locations_polar(...
%     mp.dm1.Nact/2, mp.dm1.fourier_spacing, mp.dm1.fourier_gridType, xiMin, mp.Fend.corr.Rout+1, mp.Fend.corr.ang, clocking, xiMin);


% % % % PLANNED SEARCH EFC DEFAULTS
% mp.controller = 'plannedEFC';
% mp.ctrl.dmfacVec = 1;
% mp.ctrl.log10regVec = -10:1:-3; %--log10 of the regularization exponents (often called Beta values)
% 
% %--CONTROL SCHEDULE. Columns of mp.ctrl.sched_mat are: 
%     % Column 1: # of iterations, 
%     % Column 2: log10(regularization), 
%     % Column 3: which DMs to use (12, 128, 129, or 1289) for control
%     % Column 4: flag (0 = false, 1 = true), whether to re-linearize
%     %   at that iteration.
%     % Column 5: flag (0 = false, 1 = true), whether to perform an
%     %   EFC parameter grid search to find the set giving the best
%     %   contrast .
%     % The imaginary part of the log10(regularization) in column 2 is
%     %  replaced for that iteration with the optimal log10(regularization)
%     % A row starting with [0, 0, 0, 1...] is for relinearizing only at that time
% mp.ctrl.sched_mat = [...
%     repmat([1, -3, 1, 0, 0], [4, 1]);...
%     ];
% [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);


%--Set path and filename for saved Jacobians
% mp.path.jac = % Define path to saved out Jacobians. Default if left empty is 'falco-matlab/data/jac/'
mp.jac.fn = 'jac.mat'; %'jac_iefc_test.mat'; % Name of the Jacobian file to save or that is already saved. The path to this file is set by mp.path.jac.
mp.relinItrVec = [1, 10, 20, 30]; %[];%1; %[];  %--Correction iterations at which to re-compute the Jacobian. Make an empty vector to load mp.jac.fn
% mp.relinItrVec = []; %[];%1; %[];  %--Correction iterations at which to re-compute the Jacobian. Make an empty vector to load mp.jac.fn

% Use the regular Lyot stop without pinholes for SCC
% mp.P4.compact.mask = mp.P4.compact.maskWithoutPinhole;
% mp.P4.full.mask = mp.P4.full.maskWithoutPinhole;

if ~flag_loop 
    %% Step 4: Flesh out the rest of the variables

    [mp, out] = falco_flesh_out_workspace(mp);

    %% Step 6: Perform the Wavefront Sensing and Control in FALCO

    [mp, out] = falco_wfsc_loop(mp, out);
else
%     res_arr = 7.0:0.05:7.4;
    gainfactor_arr = 0.7:0.025:1.3;
%     gainFudge_arr = 0.3:0.2:2;
%     xc_arr = [34-18.3, 34-18.3,   34-18.3, 17  ,    16.5 ];
%     yc_arr = [16,      16.5,      17,      34-18.3, 16.5 ];
    
    param_arr = gainfactor_arr;
    num_try = numel(param_arr);
    
    mp.Nitr = 3;
    InormHist_mat = zeros(num_try,mp.Nitr+1);
    betaHist_mat = zeros(num_try,mp.Nitr);
    beta_sum_mat = zeros(num_try,1);
%     InormMod_mat = zeros(num_try,mp.Nitr);
    maxStrokeDM_mat = zeros(num_try,1);

    mp0 = mp;
    for II=1:num_try
%         mp0.Fend.res = res_arr(II);
        mp0.dm1.VtoH = (4e-7*ones(mp.dm1.Nact) * 2*sqrt(2)).*gainmap*gainfactor_arr(II);%
%         mp0.est.probe.gainFudge = gainFudge_arr(II);
%         mp0.dm1.xc = xc_arr(II);              % x-center location of DM surface [actuator widths]
%         mp0.dm1.yc = yc_arr(II);               % y-center location of DM surface [actuator widths]

        if isfield(mp,'bench')
            mp = rmfield(mp,'bench');
        end
        [mp, out] = falco_flesh_out_workspace(mp0);
        mp.bench = bench;
        mp.dm1.V = zeros(mp.dm1.Nact,mp.dm1.Nact);
        [mp,out] = falco_wfsc_loop(mp,out);
        InormHist_mat(II,:) = out.InormFiberHist;
        betaHist_mat(II,:) = out.log10regHist;
%         InormMod_mat(II,:) = out.InormMod; %nm
        maxStrokeDM_mat(II) = max(max(abs(mp.dm1.V.*mp.dm1.VtoH))); %nm

        % Sum all betas
        for KK=1:mp.Nitr
            beta_sum_mat(II) =  beta_sum_mat(II) + out.log10regHist(KK);
        end
        % Plot progress
        matInorm = InormHist_mat(:,end)./InormHist_mat(:,1);
        figure(210);plot(param_arr,matInorm);title('Inorm Mat');drawnow
        figure(211);plot(param_arr,beta_sum_mat);title('BetaSum Mat');drawnow
%         figure(212);plot(res_arr,InormMod_mat(:,1)./InormHist_mat(:,1));title('InormMod Mat');drawnow
        figure(213);plot(param_arr,maxStrokeDM_mat);title('Max Stroke DM');drawnow
    end
end

return
