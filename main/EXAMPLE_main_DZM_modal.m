%% Step 1: Obtain DM command that generates dark zone
% Option 1: Generate dark zone using any method, make sure the number of
% iterations is sufficient
clear

% Step 2: Load default model parameters

EXAMPLE_defaults_try_running_FALCO


% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%--Use just 1 wavelength for initial debugging/testing of code
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

mp.Nitr = 1; %--Number of wavefront control iterations

mp.dm_ind = [1]; 
mp.est.probe.whichDM = 1; 

mp.Fend.sides = 'bottom'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);
% Option 2: Load DM command from previous experiment, load mp and out variables 

%% Step 2: Set variables for DZM

mp.Nitr = 200; % Number of iteration for EKF 
mp.dm_ind = [1]; %--DMs used in estimation/control
mp.dm_ind_static = [2]; %--DMs ONLY holding dark zone shape, not injecting drift or part of control

%%-- Variables for ekf maintenance estimation
mp.estimator = 'modal_ekf_maintenance';
mp.est.probe.Npairs = 1; 
mp.est.probe.whichDM = 1; %--Which DM is used for dither/control
mp.est.dither = 5*9.5e-3; %--std dev of dither command for random dither [V/sqtr(iter)]
mp.est.flagUseJac = true; % EKF needs the jacobian for estimation 
mp.est.read_noise = 1; %--Read noise of detector [e-]
mp.est.dark_current = 0.01; %--Dark current of detector [e-/s]
mp.est.quantum_efficiency = 0.99; %--QE of detector
mp.est.itr_ol = [1:1:mp.Nitr].'; %--"open-loop" iterations where an image is taken with initial DM command + drift command
mp.est.itr_reset = [mp.Nitr+1];
mp.est.dither_cycle_iters = 50; %--Number of unique dither commands used

%%-- Controller variables 
mp.controller = 'modal_ekf_ctrl';
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
mp.drift.magnitude = 0.01; %2*9e-4; %--std dev of random walk [V/sqrt(iter)]
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
mp.runLabel = ['modal_DZM11'];

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
