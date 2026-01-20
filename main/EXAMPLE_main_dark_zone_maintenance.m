%% Step 1: Obtain DM command that generates dark zone
% Option 1: Generate dark zone using any method, make sure the number of
% iterations is sufficient

% EXAMPLE_try_running_FALCO

% Option 2: Load DM command from previous experiment, load mp and out variables 

%% Step 2: Set variables for DZM

mp.Nitr = 3; % Number of iteration for EKF 
mp.dm_ind = [1,2]; %--DMs used in estimation/control
mp.dm_ind_static = []; %--DMs ONLY holding dark zone shape, not injecting drift or part of control

%%-- Variables for ekf maintenance estimation
mp.estimator = 'ekf_maintenance';
mp.est.probe.Npairs = 1; 
mp.est.probe.whichDM = 2; %--Which DM is used for dither/control
mp.est.dither = 9.5e-5; %--std dev of dither command for random dither [V/sqtr(iter)]
mp.est.flagUseJac = true; % EKF needs the jacobian for estimation 
mp.est.flagUseJacAlgDiff = false; % EKF needs the jacobian for estimation 
mp.est.read_noise = 1; %--Read noise of detector [e-]
mp.est.dark_current = 0.01; %--Dark current of detector [e-/s]
mp.est.itr_ol = [1:1:mp.Nitr].'; %--"open-loop" iterations where an image is taken with initial DM command + drift command
mp.est.itr_reset = [mp.Nitr+1];
mp.est.dither_cycle_iters = 50; %--Number of unique dither commands used

%%-- Controller variables 
mp.controller = 'plannedEFC';
mp.ctrl.start_iteration = 10; 
mp.ctrl.dmfacVec = 1; 
% Set EFC tikhonov parameter 
dms = str2num(strjoin(string(mp.dm_ind), '')); %--combine dm inds to be single number
mp.ctrl.sched_mat = repmat([1, -1.0, dms, 1, 0], [mp.Nitr, 1]);% 
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);
mp.relinItrVec = [1];

%%-- Drift variables 
mp.dm_drift_ind = 1;%--which DM is drifting 
mp.drift.type = 'rand_walk';%--what type of drift is happening 
mp.drift.magnitude = 9e-6; %--std dev of random walk [V/sqrt(iter)]
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
mp.runLabel = ['DZM'];

out_dir = fullfile(mp.path.config,mp.runLabel);
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
