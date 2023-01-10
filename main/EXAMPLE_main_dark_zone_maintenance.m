%% Step 1: Obtain DM command that generates dark zone
% Option 1: Generate dark zone using any method, make sure the number of
% iterations is sufficient

EXAMPLE_try_running_FALCO

% Option 2: Load DM command from previous experiment 

% mp.SeriesNum = 8;
% mp.TrialNum = 65;
% dhmStartSoln.SeriesNum = mp.SeriesNum; % Series number of previous DM solution 
% dhmStartSoln.TrialNum = mp.TrialNum; % Trial number of previous DM solution 
% dhmStartSoln.itNum = NaN; % Iteration number for previous DM solution 
% 
% load(fullfile(tb.info.OUT_DATA_DIR,['Series',num2str(dhmStartSoln.SeriesNum),'_Trial',num2str(dhmStartSoln.TrialNum)], ...
%     ['Series',num2str(dhmStartSoln.SeriesNum),'_Trial',num2str(dhmStartSoln.TrialNum),'_config.mat']))
% mp = loadPrevDMsoln(mp, dhmStartSoln, tb.info.OUT_DATA_DIR );

%% Step 2: Set variables for DZM

mp.Nitr = 100; % Number of iteration for EKF 

%%-- Variables for ekf maintenance estimation
mp.estimator = 'ekf_maintenance';
mp.est.probe.Npairs = 1; 
mp.est.probe.whichDM = 1; %--Which DM is used for control
mp.est.dither = 9.5e-4; %--std dev of dither command for random dither [V/sqtr(iter)]
mp.est.flagUseJac = true; % EKF needs the jacobian for estimation 
mp.est.read_noise = 1; %--Read noise of detector [e-]
mp.est.dark_current = 0.01; %--Dark current of detector [e-/s]
mp.est.itr_ol = [1:1:mp.Nitr].'; %--"open-loop" iterations where an image is taken with initial DM command + drift command
mp.est.itr_reset = [mp.Nitr+1];


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
mp.dm_drift_ind = mp.dm_ind;%--which DM is drifting 
mp.drift.type = 'rand_walk';%--what type of drift is happening 
mp.drift.magnitude = 9e-5; %--std dev of random walk [V/sqrt(iter)]
mp.drift.presumed_dm_std = mp.drift.magnitude; %--std dev of random walk provided to estimator, change this to account for the uncertainty of the drift magnitude

%%--
% Exposure times
tb.info.sbp_texp = 5*ones(mp.Nsbp,1); % Exposure time for each sub-bandpass (seconds)
tb.info.sbp_texp_probe = 5*ones(mp.Nsbp,1); % Exposure time for each sub-bandpass when probing (seconds)

%%-- DM settings 
mp.dm1.V_dz = mp.dm1.V; %--DM command that generates the initial dark zone
mp.dm2.V_dz = mp.dm2.V; 

mp.dm1.V_drift = zeros(mp.dm1.Nact); %--Drift injected, initialize to 0
mp.dm2.V_drift = zeros(mp.dm2.Nact);

mp.dm1.V_shift = zeros(mp.dm1.Nact); %--DM shift command for estimator reset to avoid linearization / phase wrapping errors, initialize to zero
mp.dm2.V_shift = zeros(mp.dm2.Nact);

%% Change run label and FALCO output paths
mp.runLabel = ['EKF_Series',num2str(mp.SeriesNum),'_Trial',num2str(mp.TrialNum)];

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

%% Flesh out workspace again
[mp, out] = falco_flesh_out_workspace(mp);

%% Start DZM Loop
[mp,out] = falco_wfsc_loop(mp,out);

diary off;
