%% DZM


%% Load Defaults
EXAMPLE_defaults_DST_LC_design
mp_dzm = mp;

mp_dzm.SeriesNum = 11;
mp_dzm.TrialNum = 1;

%--Use just 1 wavelength for initial debugging/testing of code
mp_dzm.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp_dzm.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp_dzm.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

mp_dzm.SeriesNum = 11;

mp_dzm.Nitr = 10; % Number of iteration for EKF 
mp_dzm.dm_ind = 1; %--DMs used in estimation/control
mp_dzm.dm_ind_static = 2; %--DMs ONLY holding dark zone shape, not injecting drift or part of control

%%-- Variables for ekf maintenance estimation
mp_dzm.estimator = 'ekf_maintenance';
mp_dzm.est.probe.Npairs = 1; 
mp_dzm.est.probe.whichDM = 1; %--Which DM is used for dither/control
mp_dzm.est.flagUseJac = true; % EKF needs the jacobian for estimation 
mp_dzm.est.read_noise = 1; %--Read noise of detector [e-]
mp_dzm.est.dark_current = 0.01; %--Dark current of detector [e-/s]
mp_dzm.est.itr_ol = [1:1:mp_dzm.Nitr].'; %--"open-loop" iterations where an image is taken with initial DM command + drift command
mp_dzm.est.itr_reset = mp_dzm.Nitr+1;

%% -- Controller variables 
mp_dzm.controller = 'plannedEFC';
% TODO: This doesn't do anything right now!
mp_dzm.ctrl.start_iteration = 10; 
mp_dzm.ctrl.dmfacVec = 1; 
% et EFC tikhonov parameter 
mp_dzm.ctrl.sched_mat = repmat([1, -1.0, 1, 1, 0], [mp_dzm.Nitr, 1]);% 
[~, mp_dzm.relinItrVec, mp_dzm.gridSearchItrVec, mp_dzm.ctrl.log10regSchedIn, mp_dzm.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp_dzm.ctrl.sched_mat);
mp_dzm.relinItrVec = 1;

%% -- Drift variables 
mp_dzm.dm_drift_ind = 1;%--which DM is drifting 
mp_dzm.drift.type = 'rand_walk';%--what type of drift is happening 

%% -- DM settings 
mp_dzm.dm1.V = out_pre.out.dm1.Vall(:, :, end); %--DM command that generates the initial dark zone
mp_dzm.dm2.V = out_pre.out.dm2.Vall(:, :, end); 
mp_dzm.dm1.V_dz = out_pre.out.dm1.Vall(:, :, end); %--DM command that generates the initial dark zone
mp_dzm.dm2.V_dz = out_pre.out.dm2.Vall(:, :, end); 

mp_dzm.dm1.V_drift = zeros(mp_dzm.dm1.Nact); %--Drift injected, initialize to 0
mp_dzm.dm2.V_drift = zeros(mp_dzm.dm2.Nact);

mp_dzm.dm1.V_shift = zeros(mp_dzm.dm1.Nact); %--DM shift command for estimator reset to avoid linearization / phase wrapping errors, initialize to zero
mp_dzm.dm2.V_shift = zeros(mp_dzm.dm2.Nact);

[mp_dzm, out_dzm] = falco_flesh_out_workspace(mp_dzm);
