%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Injects high order wavefront error drifts into the system using the DMs
% and corrects for the drifts using an extended kalman filter to estimate
% the open-loop (drifting) electric field.  The drift is corrected for
% using the electric field conjugation controller.
%
% Main parameters:
% dither - scalar, standard deviation of the dither command [V/sqrt(iter)]
% drift - scalar, standard devation of the random walk drift [V/sqrt(iter)]
% 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic parameters
%% Rename out

out_dz = out;
% need to re-initialize out?
% rename mp and re-intialize?

clear out 

%% set defaults
% TODO: need IACT defaults here
EXAMPLE_defaults_DST_LC_design

%% Step 3: Overwrite default values as desired
%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%--Use just 1 wavelength for initial debugging/testing of code
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 2;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

mp.Nitr = 10; %--Number of wavefront control iterations


%--New variables for ekf maintenance estimation:
mp.estimator = 'ekf_maintenance';
mp.est.probe.Npairs = 1;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.dither = 8e-5; %--std dev of dither command for random dither [V/sqrt(iter)]
mp.est.flagUseJac = true;

%--Controller variables
% TODO: make maint efc controller?
% mp.controller = 'maintenanceEFC';
mp.controller = 'plannedEFC';
mp.ctrl.start_iteration = 20; % controller start iteration -need to update this somehow
mp.ctrl.dmfacVec = [1];
% Set efc tikhonov parameter
mp.ctrl.sched_mat = repmat([1, -1, 1, 0, 0],[mp.Nitr,1]); % desciption explaining this can be found here: EXAMPLE_defaults_DST_LC_design.m
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%--Drift variables
mp.dm_drift_ind = [1]; %--which dms are drifting
mp.drift.type = 'rand_walk'; %--what type of drift is happening
mp.drift.magnitude = 7e-6; %--std dev of random walk [V/sqrt(iter)]

%--Initialize tb object to make things cleaner internally in sim mode
if mp.flagSim
    mp.tb.info.sbp_texp = mp.detector.tExpUnprobedVec;
    mp.tb.info.PSFpeaks = mp.detector.peakFluxVec;
end





%% Set up old commmands as V0

mp = initialize_dm_commands(mp,out_dz);

% TODO: what do I need to clear here?


%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Initialize out  (maybe not mp, need tb obj)
% flesh_out_workspace()

[mp, out] = falco_flesh_out_workspace(mp);

%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_wfsc_loop(mp, out);


%% Functions
function mp = initialize_dm_commands(mp,out)

    %--Initialize delta DM commands
    if(any(mp.dm_ind==1))
        % update to dh and V_drift and update the dm_ind for drift
        mp.dm1.V_dz = out.dm1.Vall(:,:,end); 
        
    end
    if(any(mp.dm_ind==2))
        mp.dm2.V_dz = out.dm2.Vall(:,:,end); 
       
    end

    if(any(mp.dm_drift_ind==1))
        mp.dm1.V_drift = zeros(size(mp.dm1.V)); 
    end
    if(any(mp.dm_drift_ind==2))
        mp.dm2.V_drift = zeros(size(mp.dm2.V));
    end

end

