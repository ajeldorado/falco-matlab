clear
clc

% Load and add current folder

current_folder = pwd;
addpath(genpath(current_folder));
% addpath(genpath("C:\Users\skmk9\Downloads"));

%% Load dm_commands
out_pre_config = load("C:\Users\skmk9\Downloads\Series0001_Trial0001_LC_simple_2DM48_z1_IWA2.8_OWA10_5lams550nm_BW10_plannedEFC_config.mat");
out_pre = load("C:\Users\skmk9\Downloads\Series0001_Trial0001_LC_simple_2DM48_z1_IWA2.8_OWA10_5lams550nm_BW10_plannedEFC_snippet.mat");

init_command_dm1 = out_pre.out.dm1.Vall(:,:,end);  % 2D voltage map for dm1 at iteration Itr
init_command_dm2 = out_pre.out.dm2.Vall(:,:,end);  % 2D voltage map for dm2 at iteration Itr
% Load the second to last command
delta_dm1 = abs(init_command_dm1-out_pre.out.dm1.Vall(:,:,end-1));
delta_dm2 = abs(init_command_dm2-out_pre.out.dm2.Vall(:,:,end-1));
dm1.V_dz = init_command_dm1;
dm2.V_dz = init_command_dm2;

% delta_dm_mean = mean([std(delta_dm1(init_command_dm1 ~= 0)), std(delta_dm2(init_command_dm2 ~= 0))]);
% if size(dm_ind,2) > 1
%     mp.delta_dm_mean = mean([mean(mp.delta_dm1(mp.init_command_dm1 ~= 0)), mean(mp.delta_dm2(mp.init_command_dm2 ~= 0))]);

% elseif any(mp.dm_ind == 1)

delta_dm_mean = mean([mean(delta_dm1(init_command_dm1 ~= 0))]);

% elseif any(mp.dm_ind == 2)
% 
%     mp.delta_dm_mean = mean([mean(mp.delta_dm2(mp.init_command_dm2 ~= 0))]);
% 
% end

DZMExperimentTemplate
PWPExperimentTemplate

drifts = [delta_dm_mean/10];

mp_dzm_arr = {};
out_dzm_arr = {};
mp_arr = {};
out_arr = {};

base_path = mp_dzm.path.config;

%% Step 2: Set variables for DZM
for i = 1:length(drifts)
    mp_dzm.drift.magnitude = drifts(i); %--std dev of random walk [V/sqrt(iter)]
    mp_dzm.drift.presumed_dm_std = mp_dzm.drift.magnitude; %--std dev of random walk provided to estimator, change this to account for the uncertainty of the drift magnitude
    dither_vec = [5*drifts(i), 10*drifts(i)];
    for ii = dither_vec
        mp_dzm.est.dither = ii;
        
        name = ['diary_', num2str(mp_dzm.TrialNum),'.txt'];

        mp_dzm.runLabel = ['DZM_', 'Series',num2str(mp_dzm.SeriesNum,'%04d'),'_Trial',num2str(mp_dzm.TrialNum,'%04d_')];

        out_dir = fullfile(base_path, mp_dzm.runLabel);
        mp_dzm.path.config = out_dir;  %--Location of *config.mat and *snippet.mat output files
        mp_dzm.path.ws = out_dir; % Location for (mostly) complete workspace from end of trial, *_all.mat
        
        % Directory to save data
        if(~exist(out_dir, 'dir'))
            mkdir(out_dir);
        else
            disp('Output directory already exists. Press enter to continue and overwrite contents.');
            pause
        end
        mp_dzm.diaryfile = fullfile(out_dir,name);
        
        tb = TestbedExperimentTemplate(mp_dzm);

        diary(mp_dzm.diaryfile)
        
        %% Start DZM Loop
        [mp_dzm_out, out_dzm_out] = falco_wfsc_loop(mp_dzm,out_dzm);
        
        mp_dzm_arr{end+1} = mp_dzm_out;
        out_dzm_arr{end+1} = out_dzm_out;
        diary off;
        
        %% Update Trial Numbers
        mp.TrialNum = mp_dzm.TrialNum+1;
        mp_dzm.TrialNum = mp.TrialNum+1;

        name = ['diary_', num2str(mp.TrialNum),'.txt'];
        mp.diaryfile = fullfile(out_dir,name);
        diary(mp.diaryfile)
        %% Start PWP Loop
        tb = TestbedExperimentTemplate(mp);
        % ev = resetPins(ev, mp_dzm, out_pre);
        [mp_out, out_pwp] = falco_wfsc_loop(mp,out);

        mp_arr{end+1} = mp_out;
        out_arr{end+1} = out_pwp;
        diary off;

        % ev = resetPins(ev, mp, out_pre);

    end
end

plotParamScanData(mp_dzm_arr);


function ev = resetPins(ev, mp, out_pre)
    % Initialize pinned actuator check
    ev.dm1.initial_pinned_actuators = out_pre.out.dm1.pinned{end};
    if any(mp.dm_ind == 2); ev.dm2.initial_pinned_actuators = out_pre.out.dm2.pinned{end}; end
    
    ev.dm1.new_pinned_actuators = [];
    ev.dm2.new_pinned_actuators = [];
    ev.dm1.act_ele_pinned = [];
    ev.dm2.act_ele_pinned = [];
    
    ev.exit_flag = false;
end