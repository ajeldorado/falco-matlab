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
flag_smallDH4smf = false; % Arielle!
flag_efc360is = false;
mp.flag_timeMaya = false;
mp.flag_adjusttime = true;
%% Step 2: Load default model parameters
flag_svc = false;
flag_svc_dimple = false;

if mp.flag_timeMaya
    Time_total_start = tic(); %start a tic/toc for total runnning time
end 

if flag_svc
    lambda0 = 760*1e-9;
elseif flag_svc_dimple
    lambda0 = 760*1e-9;
else
    lambda0 = 780*1e-9;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsbp = 1; % Number of sub-bandpasses, AKA number of wavelengths in full bandpass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NDfilter_FWpos = 2;
FWpos = 1;  

if FWpos~=NDfilter_FWpos 
    NDfilter_cal = 28.6082 * 4.9367 * 4.4753 * 4.6956;%23.16; %>1
%     NDfilter_cal = 28.6082 * 4.9367 * 4.4753;% * 4.6956;%23.16; %>1
else
    NDfilter_cal = 1;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tint_offAxis = 5e-2*ones(Nsbp,1); %; % 1e-3; %Integration time for off-axis PSF
mp.tint_efc = 1e-2*ones(Nsbp,1); %2.5e-0*ones(Nsbp,1); %30; %1e-1; %30; %35e-0;%1e-1;%[2,1,1,0.5,1]*0.05; % for the control and evalution
mp.tint_est = 1e-2*ones(Nsbp,1);%5e-1; % time of integration for the estimation/sensing
% tint_offAxis = 5e-2*ones(Nsbp,1); %; % 1e-3; %Integration time for off-axis PSF
% mp.tint_efc = 2.5e-0*ones(Nsbp,1); %30; %1e-1; %30; %35e-0;%1e-1;%[2,1,1,0.5,1]*0.05; % for the control and evalution
% mp.tint_est = 2.5e-1*ones(Nsbp,1);%5e-1; % time of integration for the estimation/sensing

% tint_offAxis(1) = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Source parameters
bench.info.source = 'nkt';%'laser';% 

mp.Nsbp = Nsbp;
if Nsbp>1
    mp.fracBW = 0.1;  % bandwidth of experiment
    mp.sbp_centers = lambda0*linspace( 1-mp.fracBW/2,1+mp.fracBW/2,Nsbp);
    bench.info.sbp_width = 8e-9*ones(Nsbp,1);%[14,14,14,14,14]*1e-9;% 3e-9;%-  %--Width of each sub-bandpass on testbed (meters)
%     bench.info.sbp_width([4:6]) = 6e-9; %
else
    mp.fracBW = 0.013; % bandwidth of experiment; when Nsbp is 1 --> narrowband or monochromatic  
%     mp.fracBW = 0.1;
    mp.sbp_centers = lambda0;
    bench.info.sbp_width = [8]*1e-9; % bandpass of the sub-bandpass [nm]
%     bench.info.sbp_width = [78]*1e-9; % bandpass of the sub-bandpass [nm]
    mp.est.probe.gainFudge = 1;
    if flag_efc360is
        mp.est.probe.gainFudge = 3;
    end
end

%%
mp.flagFiber = false;
flag_apodizer_avc = false;
EXAMPLE_defaults_HCST_VVC % Re-use SCC example's configuration
mp.flagSim = false;      %--Simulation or not
mp.testbed = 'HCST';

% mp.Nitr = 49;
%% HCST preparation
% Flags; all flags false --> regular EFC run
flagFieldStop = false; % Arielle!

% path2flatMap = [bench.info.HCST_DATA_DIR,'zygo_flat/2019Jun26/'];
% path2flatMap = [bench.info.HCST_DATA_DIR,'is_results/2023Sep15/'];
% path2flatMap =[bench.info.HCST_DATA_DIR,'pr_results/2023Jun13/'];
% path2flatMap = [bench.info.HCST_DATA_DIR,'is_results/2023Dec02/'];
% path2flatMap = [bench.info.HCST_DATA_DIR,'pr_results/2023Nov29/'];
% path2flatMap = [bench.info.HCST_DATA_DIR,'pr_results/2024Apr09/'];
path2flatMap = [bench.info.HCST_DATA_DIR,'pr_results/2024Oct04/'];
% path2flatMap = [bench.info.HCST_DATA_DIR,'pr_results/2024Apr12/'];

% path2startingMap = [bench.info.HCST_DATA_DIR,'zygo_flat/2019Jun26/'];
% path2startingMap = [bench.info.HCST_DATA_DIR,'pr_results/2023Nov29/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'is_results/2023Sep15/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'pr_results/2023Jun13/'];
% path2startingMap = [bench.info.HCST_DATA_DIR,'pr_results/2024Apr09/'];
path2startingMap = [bench.info.HCST_DATA_DIR,'pr_results/2024Oct04/'];

% path2startingMap =[bench.info.HCST_DATA_DIR,'efc_results/2023Aug14_sol360/'];
%  path2startingMap =[bench.info.HCST_DATA_DIR,'efc_results/2023Aug18_sol/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'efc_results/2023Jun29_sol/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'pr_results/2023Jun13/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'efc_results/2024Aug13/'];
% path2startingMap = [bench.info.HCST_DATA_DIR,'is_results/2023Dec02/'];
% path2startingMap = [bench.info.HCST_DATA_DIR,'pr_results/2023Nov29/'];

frameSize = 180;

if flag_efc360is
    mp.Nitr = 9;
    mp.est.probe.gainFudge = 0.2;
    frameSize = 180;
    Nsbp = 1;
    mp.tint_efc = 5e-1*ones(Nsbp,1); %30; %1e-1; %30; %35e-0;%1e-1;%[2,1,1,0.5,1]*0.05; % for the control and evalution
    mp.tint_est = 2.5e-1*ones(Nsbp,1);%5e-1; % time of integration for the estimation/sensing
    flagFieldStop = false;
end

datelabel = datestr(now,'yyyymmmddTHHMM');
label_progress = [num2str(Nsbp),'lams',num2str(round(1e9*lambda0)),'nm_BW',num2str(mp.fracBW*100),'_VVC_',datelabel];

falcoPreparev2

mp.bench = bench;
% bench.info.sbp_width = sbp_width ; %--Width of each sub-bandpass on testbed (meters)
bench.info.sbp_texp = mp.tint_efc;% Exposure time for each sub-bandpass (seconds)
bench.info.PSFpeaks = PSFpeaks;% counts per second 


%% Step 3: Overwrite default values as desired

%--Whether to perform a model-based (instead of empirical) grid search for the controller
mp.ctrl.flagUseModel = false; 

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = false; %% MAYA!


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
mp.relinItrVec = [1]; %[];%1; %[];  %--Correction iterations at which to re-compute the Jacobian. Make an empty vector to load mp.jac.fn
% mp.relinItrVec = []; %[];%1; %[];  %--Correction iterations at which to re-compute the Jacobian. Make an empty vector to load mp.jac.fn

% Use the regular Lyot stop without pinholes for SCC
% mp.P4.compact.mask = mp.P4.compact.maskWithoutPinhole;
% mp.P4.full.mask = mp.P4.full.maskWithoutPinhole;

if ~flag_loop 
    %% Step 4: Flesh out the rest of the variables

    [mp, out] = falco_flesh_out_workspace(mp);

    %% Step 6: Perform the Wavefront Sensing and Control in FALCO

    [mp, out] = falco_wfsc_loop(mp, out);
%     hcst_disconnectDevices(bench, true, false);
else
%      zrot_arr = 0.1:0.05:0.3;
%     res_arr = 7.05:0.01:7.1;
%     gainfactor_arr = 0.4:0.025:0.6;
%     gainFudge_arr = 0.3:0.2:2;
%     xc_arr = [16.,    16.,         16.,         16.+0.5, 16.+0.5,    16.+0.5,    16.-0.5, 16.-0.5,     16.-0.5];
%     yc_arr = [16., 16.+0.5, 16.-0.5, 16.,   16.+0.5, 16.-0.5, 16.,    16.+0.5, 16.-0.5];
%     xc_arr = [16.25,  16.25, 16.25,  16.25, 16.25, 16.25, 16.25, 16.25];
%     yc_arr = [15.,  15.5,   16.,  16.5,   17 , 17.5,  18, 18.5];
%     xc_arr = [16.25,  16.25, 16.25,  16.25, 16.25, 16.25, 16.25, 16.25];
%     yc_arr = [15.,  15.5,   16.,  16.5,   17 , 17.5,  18, 18.5];
    xc_arr = repmat([16, 16.5, 17., 17.5], 1, 4);
    yc_arr = repelem([16, 16.5, 17, 17.5], 4) ;

%     param_arr = zrot_arr;
    param_arr = 1:numel(xc_arr);
    num_try = numel(param_arr);
    
    mp.Nitr = 5;
    InormHist_mat = zeros(num_try,mp.Nitr+1);
    betaHist_mat = zeros(num_try,mp.Nitr);
    beta_sum_mat = zeros(num_try,1);
%     InormMod_mat = zeros(num_try,mp.Nitr);
    maxStrokeDM_mat = zeros(num_try,1);

    mp0 = mp;
    for II=1:num_try
%         mp0.Fend.res = res_arr(II);%*mp.lambda0/(785e-9);
%         mp0.dm1.VtoH = (4e-7*ones(mp.dm1.Nact) * 2*sqrt(2)).*gainmap*gainfactor_arr(II);%
%         mp0.est.probe.gainFudge = gainFudge_arr(II);
        mp0.dm1.xc = xc_arr(II);              % x-center location of DM surface [actuator widths]
        mp0.dm1.yc = yc_arr(II);               % y-center location of DM surface [actuator widths]
%         mp.dm1.zrot = zrot_arr(II);
        
        if isfield(mp,'bench')
            mp = rmfield(mp,'bench');
        end
        mp0.Fend.FOV = mp.Fend.Nxi/2/mp0.Fend.res; 
        [mp, out] = falco_flesh_out_workspace(mp0);
        mp.bench = bench;
        mp.dm1.V = zeros(mp.dm1.Nact,mp.dm1.Nact);
        [mp,out] = falco_wfsc_loop(mp,out);
        InormHist_mat(II,:) = out.InormHist;
        betaHist_mat(II,:) = out.log10regHist;
%         InormMod_mat(II,:) = out.InormMod; %nm
        maxStrokeDM_mat(II) = max(max(abs(mp.dm1.V.*mp.dm1.VtoH))); %nm

        % Sum all betas
        for KK=1:mp.Nitr
            beta_sum_mat(II) =  beta_sum_mat(II) + out.log10regHist(KK);
        end
        % Plot progress
        matInorm = InormHist_mat(:,end); %./InormHist_mat(:,1);
        figure(210);plot(param_arr,matInorm);title('Inorm Mat');drawnow
        figure(211);plot(param_arr,beta_sum_mat);title('BetaSum Mat');drawnow
%         figure(212);plot(res_arr,InormMod_mat(:,1)./InormHist_mat(:,1));title('InormMod Mat');drawnow
        figure(213);plot(param_arr,maxStrokeDM_mat);title('Max Stroke DM');drawnow
    end
end


if mp.flag_timeMaya
    Time_total = toc(Time_total_start); %end of tic/toc for total running time
    fprintf('Total Running time for %d iterations: %.2f seconds\n', Itr, Time_total);
end 

return
%%



