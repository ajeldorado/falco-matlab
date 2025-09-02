%% Step 1: Obtain DM command that generates dark zone
% Option 1: Generate dark zone using any method, make sure the number of
% iterations is sufficient

EXAMPLE_try_running_FALCO

% Option 2: Load DM command from previous experiment, load mp and out variables 

%% Step 2: Set variables for DZM

mp.Nitr = 100; % Number of iteration for EKF 
mp.dm_ind = [1,2]; %--DMs used in estimation/control
mp.dm_ind_static = []; %--DMs ONLY holding dark zone shape, not injecting drift or part of control

%%-- Variables for ekf maintenance estimation
mp.estimator = 'aekf_maintenance';
mp.est.probe.Npairs = 1; 
mp.est.probe.whichDM = 2; %--Which DM is used for dither/control
mp.est.dither = 9.5e-5; %--std dev of dither command for random dither [V/sqtr(iter)]
mp.est.flagUseJac = true; % EKF needs the jacobian for estimation 
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

%%-- Planet injection parameters - just the mp.star setup
planet_inds = [2];
inj_pos = [1,5];
mp.star.xiOffsetVec = [0, inj_pos(1)];
mp.star.etaOffsetVec = [0, inj_pos(2)];
mp.star.weights = [1,1e-7];
num_rolls = 4;
roll_count = 1;
roll_angle_deg = 25;
roll_angle = roll_angle_deg * pi / 180;
current_roll_angle = 0;
mp.star.sep_rad = sqrt((mp.star.etaOffsetVec ).^2 + (mp.star.xiOffsetVec ).^2);
mp.star.pos_angle = atan2( mp.star.etaOffsetVec ,  mp.star.xiOffsetVec );

mp.compact.star.count = length(mp.star.weights);
mp.compact.star.xiOffsetVec = mp.star.xiOffsetVec;
mp.compact.star.etaOffsetVec = mp.star.etaOffsetVec;
mp.compact.star.weights = mp.star.weights;

%% CRITICAL FIX: Set up Jacobian star weights BEFORE fleshing out workspace
% This must match the number of stars to prevent indexing errors
mp.jac.star.weights = mp.star.weights; % Should be [1, 1e-7] for 2 stars

%% Enable incoherent estimation in EKF
% Set flag to use 3-element state (real, imag, incoherent) instead of 2-element (real, imag)
mp.est.flagIncoherent = true; % Enable incoherent intensity estimation

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
mp.runLabel = ['DZM_with_planet_aekf'];

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

%% Debug check before fleshing out workspace
fprintf('Before falco_flesh_out_workspace:\n');
fprintf('  mp.star.count: %d\n', length(mp.star.weights));
fprintf('  mp.compact.star.count: %d\n', mp.compact.star.count);
fprintf('  mp.jac.star.weights length: %d\n', length(mp.jac.star.weights));
fprintf('  mp.star.weights: [%s]\n', num2str(mp.star.weights));
fprintf('  mp.jac.star.weights: [%s]\n', num2str(mp.jac.star.weights));
fprintf('  Incoherent estimation enabled: %s\n', mat2str(mp.est.flagIncoherent));

%% Flesh out workspace again
[mp, out] = falco_flesh_out_workspace(mp);

%% Make sure spatial weighting has one column per star used in the Jacobian.
% Many configs set WspatialVec as a single column (1 star). With multi-star
% (e.g., star + injected planet), falco_apply_spatial_weighting_to_Jacobian
% will index (:, iStar), so we must replicate columns.
if isfield(mp,'WspatialVec') && ~isempty(mp.WspatialVec)
    if size(mp.WspatialVec,2) == 1 && mp.compact.star.count > 1
        mp.WspatialVec = repmat(mp.WspatialVec, 1, mp.compact.star.count);
        fprintf('Replicated WspatialVec to size [%d x %d] for %d stars.\n', ...
            size(mp.WspatialVec,1), size(mp.WspatialVec,2), mp.compact.star.count);
    end
else
    % If WspatialVec wasn't created, build a flat weighting over correction region.
    % Use correction mask length; make one column per star.
    if isfield(mp,'Fend') && isfield(mp.Fend,'corr') && isfield(mp.Fend.corr,'mask')
        baseW = double(mp.Fend.corr.mask(:));
        mp.WspatialVec = repmat(baseW, 1, mp.compact.star.count);
        fprintf('Created WspatialVec of size [%d x %d] for %d stars.\n', ...
            size(mp.WspatialVec,1), size(mp.WspatialVec,2), mp.compact.star.count);
    end
end

%% Set up mp.star.count for multi-star simulation
mp.star.count = length(mp.star.weights); % Should be 2 (star + planet)

%% Final debug check
fprintf('After setup:\n');
fprintf('  mp.star.count: %d\n', mp.star.count);
fprintf('  mp.compact.star.count: %d\n', mp.compact.star.count);
if isfield(mp, 'WspatialVec')
    fprintf('  WspatialVec size: [%d x %d]\n', size(mp.WspatialVec,1), size(mp.WspatialVec,2));
end

%% Test the multi-star simulation
disp('Testing multi-star image generation...')
I_test = falco_get_sim_sbp_image(mp, 1);
disp(['Image generated with ', num2str(mp.star.count), ' stars'])
figure; imagesc(log10(I_test)); colorbar; title('Multi-star test image'); axis equal; axis tight;

%% Start DZM Loop
[mp,out] = falco_wfsc_loop(mp,out);

diary off;