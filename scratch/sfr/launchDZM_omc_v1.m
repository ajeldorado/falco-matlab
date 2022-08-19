% Script to launch dark hole maintainance after running falco_wfsc_loop the
% traditional way 
falco_omc_paths;
tb.info.OUT_DATA_DIR = data_dir;
%% Use a previous DM solution 
% mp.SeriesNum = 8;
% mp.TrialNum = 65;
dhmStartSoln.SeriesNum = mp.SeriesNum; % Series number of previous DM solution 
dhmStartSoln.TrialNum = mp.TrialNum; % Trial number of previous DM solution 
dhmStartSoln.itNum = NaN; % Iteration number for previous DM solution 

load(fullfile(tb.info.OUT_DATA_DIR,['Series',num2str(dhmStartSoln.SeriesNum),'_Trial',num2str(dhmStartSoln.TrialNum)], ...
    ['Series',num2str(dhmStartSoln.SeriesNum),'_Trial',num2str(dhmStartSoln.TrialNum),'_config.mat']))
mp = loadPrevDMsoln(mp, dhmStartSoln, tb.info.OUT_DATA_DIR );
%mp = loadPrevEsens(mp, dhmStartSoln, tb.info.OUT_DATA_DIR );

%% 

mp.Nitr = 450; % Number of iteration for EKF 

%%-- New variables for ekf maintenance estimation
mp.estimator = 'ekf_maintenance';
mp.est.probe.Npairs = 1; 
mp.est.probe.whichDM = 1; 
mp.est.dither = 0.025; %--std dev of dither command for random dither [V/sqtr(iter)]
mp.est.flagUseJac = true; % EKF needs the jacobian for estimation 
mp.est.read_noise = 1; % e-
mp.est.dark_current = 0.01; %e-/s
mp.est.itr_ol = [1:200:mp.Nitr,mp.Nitr].'; % "open-loop" iterations 
mp.est.itr_reset = [500,600];


%%-- Controller variables 
mp.controller = 'plannedEFC';
% TODO: This doesn't do anything right now!
mp.ctrl.start_iteration = 10; 
mp.ctrl.dmfacVec = 1; 
% Set EFC tikhonov parameter 
mp.ctrl.sched_mat = repmat([1, -1.0, 1, 1, 0], [mp.Nitr, 1]);% 
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);
mp.relinItrVec = [];

%%-- Drift variables 
mp.dm_drift_ind = 1;%--which DM is drifting 
mp.drift.type = 'rand_walk';%--what type of drift is happening 
mp.drift.magnitude = 0.0025; %--std dev of random walk [V/sqrt(iter)]
mp.drift.presumed_dm_std =0.0025; 

%%--
% Exposure times
tb.info.sbp_texp = 5*ones(mp.Nsbp,1); % Exposure time for each sub-bandpass (seconds)
tb.info.sbp_texp_probe = 5*ones(mp.Nsbp,1); % Exposure time for each sub-bandpass when probing (seconds)

%%-- DM settings 
mp.dm1.V_dz = mp.dm1.V; 
mp.dm2.V_dz = mp.dm2.V; 

mp.dm1.V_drift = zeros(mp.dm1.Nact);
mp.dm2.V_drift = zeros(mp.dm2.Nact);

mp.dm1.V_shift = zeros(mp.dm1.Nact);
mp.dm2.V_shift = zeros(mp.dm2.Nact);

%% Change run label and FALCO output paths
% mp.TrialNum = 66;
mp.runLabel = ['EKF_Series_omc_',num2str(mp.SeriesNum),'_Trial',num2str(mp.TrialNum)];

out_dir = fullfile(tb.info.OUT_DATA_DIR,mp.runLabel);
mp.path.config = out_dir;  %--Location of *config.mat and *snippet.mat output files
mp.path.ws = out_dir; % Location for (mostly) complete workspace from end of trial, *_all.mat

% Directory to save data
if(~exist(out_dir, 'dir'))
    mkdir(out_dir);
else
    disp('Output directory already exists. Press enter to continue and overwrite contents.');
    pause
end

mp.diaryfile = fullfile(out_dir,'diary.txt');
diary(mp.diaryfile)

%%

[mp, out] = falco_flesh_out_workspace(mp);
mp.tb = tb; % This hides all of the testbed parameters in the model params 
            % so custom functions internal to Falco can use it.
if mp.flagSim
    mp.tb.info.sbp_texp = mp.detector.tExpUnprobedVec;
    mp.tb.info.PSFpeaks = mp.detector.peakFluxVec;
end

%%
[mp,out] = falco_wfsc_loop(mp,out);

%%

diary off;
