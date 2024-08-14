%% Load files
basepath = "";
SeriesNum = 001;
TrialNum = 001;

load(fullfile(basepath, ['Series',num2str(SeriesNum),'_Trial',num2str(TrialNum),'_config.mat']))

% Testbed
if mp.testbed
    StartSoln.SeriesNum = mp.SeriesNum; % Series number of previous DM solution
    StartSoln.TrialNum = mp.TrialNum; % Trial number of previous DM solution

    % Get the previous of the previous DM solution
    StartSoln.ItNum = mp.Nitr-1; % Iteration number of previous of the previous DM solution
    mp = loadPrevDMsoln(mp, StartSoln, tb.info.OUT_DATA_DIR );

    out.dm1.Vall = zeros(mp.dm1.Nact, mp.dm1.Nact, 2);
    out.dm2.Vall = zeros(mp.dm2.Nact, mp.dm2.Nact, 2);

    out.dm1.Vall(:, :, 1) = mp.dm1.V;
    out.dm2.Vall(:, :, 1) = mp.dm2.V;

    StartSoln.itNum = NaN; % Iteration number for previous DM solution
    mp = loadPrevDMsoln(mp, StartSoln, tb.info.OUT_DATA_DIR );

    out.dm1.Vall(:, :, 2) = mp.dm1.V;
    out.dm2.Vall(:, :, 2) = mp.dm2.V;

% Simulation
else
    load(fullfile(basepath, ['Series',num2str(SeriesNum),'_Trial',num2str(TrialNum),'_all.mat']))
end

%% Initialize pinned act list
% Initialize pinned actuator check
ev.dm1.initial_pinned_actuators = mp.dm1.pinned;
if any(mp.dm_ind == 2); ev.dm2.initial_pinned_actuators = mp.dm2.pinned; end

ev.dm1.new_pinned_actuators = [];
ev.dm2.new_pinned_actuators = [];
ev.dm1.act_ele_pinned = [];
ev.dm2.act_ele_pinned = [];

ev.exit_flag = false;

%% Select the specific iteration index
mp.Itr = mp.Nitr; 
mp.dm_ind = [1];
mp.dm_ind_static = [2];

%% Setting initial DM commands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Access the 2D voltage maps for specified iteration
mp.init_command_dm1 = out.dm1.Vall(:,:,end);  % 2D voltage map for dm1 at iteration Itr
mp.init_command_dm2 = out.dm2.Vall(:,:,end);  % 2D voltage map for dm2 at iteration Itr

% Load the second to last command
mp.delta_dm1 = abs(mp.init_command_dm1-out.dm1.Vall(:,:,end-1));
mp.delta_dm2 = abs(mp.init_command_dm2-out.dm2.Vall(:,:,end-1));

mp.dm1.V_dz = mp.init_command_dm1;
mp.dm2.V_dz = mp.init_command_dm2;

% delta_dm_mean = mean([std(delta_dm1(init_command_dm1 ~= 0)), std(delta_dm2(init_command_dm2 ~= 0))]);
if size(mp.dm_ind,2) > 1
    out.delta_dm_mean = mean([mean(mp.delta_dm1(mp.init_command_dm1 ~= 0)), mean(mp.delta_dm2(mp.init_command_dm2 ~= 0))]);

elseif any(mp.dm_ind == 1)

    out.delta_dm_mean = mean([mean(mp.delta_dm1(mp.init_command_dm1 ~= 0))]);

elseif any(mp.dm_ind == 2)

    out.delta_dm_mean = mean([mean(mp.delta_dm2(mp.init_command_dm2 ~= 0))]);

end

fprintf('Successfully loaded initial DM commands.\n');

%% Setting initial variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.est.dithers = out.delta_dm_mean*[0.1, 1, 3, 5, 10, 15, 20, 30, 50]; % ideally 3 dithers
% out.est.dithers = [1e-10, 1e-5, 1e-3, 1];
out.num_dithers = numel(out.est.dithers); % set to the same index number of out.est.dithers
out.num_iterations = 25;
% num_subbands = 1;
out.iSubband = ceil(mp.Nsbp / 2); % Pick the middle wavelength
fprintf('Now setting necessary initial variables.\n');

%% Set total command for estimator image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: need to save these commands for each iteration separately
fprintf('Now entering first for loop. (give statistics of each dither)\n');


%% Create an array for mp.output arrays
% contrasts --> # dithers (rows) vs # iterations (columns)

% this needs to be 2 dimensional
out.contrasts = zeros(out.num_dithers, out.num_iterations);
out.dithers_nm_1 = zeros(out.num_dithers, out.num_iterations);
if any(mp.dm_ind == 2)
out.dithers_nm_2 = zeros(out.num_dithers, out.num_iterations);
end

%% For loops %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One loop that goes over the iterations
for i = 1:out.num_dithers
    dither = out.est.dithers(i); % assuming only one dither value, otherwise use i if there are multiple dithers
    fprintf(['on dither ',num2str(i)])
    % Second loop that will give statistics of each dither
    for j = 1:out.num_iterations
        fprintf('Now inside second for loop.\n');
        
        V1 = zeros(mp.dm1.NactTotal,1);
        V2 = zeros(mp.dm2.NactTotal,1);

        % Calculate the dither commands
        if any(mp.dm_ind == 1)
            V1(mp.dm1.act_ele) = normrnd(0,out.est.dithers(i),size(mp.dm1.act_ele));
            DM1Vdither = reshape(V1,[mp.dm1.Nact,mp.dm1.Nact]);
        else
            DM1Vdither = zeros(size(mp.dm1.V));
        end

        if any(mp.dm_ind == 2)
            V2(mp.dm2.act_ele) = normrnd(0,out.est.dithers(i),size(mp.dm2.act_ele));
            DM2Vdither = reshape(V2,[mp.dm2.Nact,mp.dm2.Nact]);
        else
            DM2Vdither = zeros(size(mp.dm2.V));
        end

        % Generate command to apply to DMs
        if any(mp.dm_ind == 1)
            % Note falco_set_constrained_voltage does not apply the command to the DM
            mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + DM1Vdither); 
        elseif any(mp.dm_ind_static == 1)
            mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz);
        end

        if any(mp.dm_ind == 2)
            mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + DM2Vdither); 
        elseif any(mp.dm_ind_static == 2)
            mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz);
        end

        % Do safety check to make sure no actuators are pinned
        % ev = pinned_act_safety_check(mp, ev);
        % closed_loop_command = dither + efc_command + get_dm_command_vector(mp, mp.dm1.V_shift, mp.dm2.V_shift);

        % Debugging statements before y_measured assignment 
        disp('Assigning y_measured:');
        disp(['mp.Fend.corr.Npix: ', num2str(mp.Fend.corr.Npix)]);
        disp(['mp.Nsbp: ', num2str(mp.Nsbp)]);

        % Assign y_measured
        y_measured = zeros(mp.Fend.corr.Npix, mp.Nsbp);
        disp('y_measured assigned successfully');


        [mp, ev] = pinned_act_safety_check(mp,ev);
        % Only take image if no actuators in act_ele will get pinned
        if ~ev.exit_flag
            % Call falco_get_sbp_image for each subband
            disp('Calling falco_get_sbp_image...');
            I = falco_get_sbp_image(mp, out.iSubband); % iSubband : index of subband for which to take the image
            disp('falco_get_sbp_image called successfully');

            % Need to convert if it's not in contrast units
            % I0 = ev.imageArray / peak to put it into contrast units % ev is the estimator variable... (internal)
            mean_contrast = mean(I(mp.Fend.score.mask)); % this is the dark zone mask
            disp(['Mean contrast for subband ', num2str(out.iSubband), ': ', num2str(mean_contrast)]);
            
            % Generate command to apply to DMs
            if any(mp.dm_ind == 1)
                % Note falco_set_constrained_voltage does not apply the command to the DM
                mp.dm1 = falco_set_constrained_voltage(mp.dm1, DM1Vdither); 
            end
    
            if any(mp.dm_ind == 2)
                mp.dm2 = falco_set_constrained_voltage(mp.dm2, DM2Vdither); 
            end
            
            
            %% STILL NEED TO DO:
            if any(mp.dm_ind == 1)
            
                DM1Vdither_nm = falco_calc_act_height_from_voltage(mp.dm1)*10^9; % convert to nano meters
                % Calculate standard deviation of the dithered commands (excluding zeros)
                out.dithers_nm_1(i, j) = std(DM1Vdither_nm(mp.init_command_dm1 ~= 0));
            end
            if any(mp.dm_ind == 2)
                DM2Vdither_nm = falco_calc_act_height_from_voltage(mp.dm2)*10^9;
                out.dithers_nm_2(i, j) = std(DM2Vdither_nm(mp.init_command_dm2 ~= 0));
            end
            % Store results
            out.contrasts(i, j) = mean_contrast;
        end
        
    end
end

%% Plotting
% four things:
%   1. contrast per iteration for each of the dithers
figure;
for i = 1:out.num_dithers
    semilogy(1:out.num_iterations, out.contrasts(i, :), '-o');
    hold on;
end
title('Contrast per iteration for each of the dithers');
xlabel('Iteration');
ylabel('Contrast');
legend(arrayfun(@(x) ['Dither ' num2str(x)], 1:out.num_dithers, 'UniformOutput', false));
hold off;

%   2. dither std nms vs. iteration for each of the dithers
figure;
for i = 1:out.num_dithers
    plot(1:out.num_iterations, out.dithers_nm_1(i, :), '-o');
    hold on;
end
title('DM1 Dither std nms vs. iteration for each of the dithers');
xlabel('Iteration');
ylabel('Standard Deviation (nm)');
legend(arrayfun(@(x) ['Dither ' num2str(x)], 1:out.num_dithers, 'UniformOutput', false));
hold off;

if any(mp.dm_ind == 2)
figure;
for i = 1:out.num_dithers
    plot(1:out.num_iterations, out.dithers_nm_2(i, :), '-o');
    hold on;
end
title('DM2 Dither std nms vs. iteration for each of the dithers');
xlabel('Iteration');
ylabel('Standard Deviation (nm)');
legend(arrayfun(@(x) ['Dither ' num2str(x)], 1:out.num_dithers, 'UniformOutput', false));
hold off;
end

%   3. take the averages along the dimensions of both the contrasts and the dither std
%                   nms ver the iteration axis, then you plot the mean std
%                   dithers (as x-axis) and mean contrast (as y-axis)
mean_contrasts = mean(out.contrasts, 2);
mean_dithers_nm_1 = mean(out.dithers_nm_1, 2);


figure;
plot(mean_dithers_nm_1, mean_contrasts, '-o');
title('Mean Contrast vs. Mean DM1 Dither std nms');
xlabel('Mean DM1 Dither std (nm)');
xlim([min(mean_dithers_nm_1)/2, max(mean_dithers_nm_1)*1.5]);
ylabel('Mean Contrast');
grid on;

if any(mp.dm_ind == 2)
mean_dithers_nm_2 = mean(dithers_nm_2, 2);
figure;
plot(mean_dithers_nm_2, mean_contrasts, '-o');
title('Mean Contrast vs. Mean DM2 Dither std nms');
xlabel('Mean DM2 Dither std (nm)');
xlim([min(mean_dithers_nm_2)/2, max(mean_dithers_nm_2)*1.5]);
ylabel('Mean Contrast');
grid on;
end

%   4. mean of the mean contrasts (y-axis) vs. dithers in voltages (x-axis)
mean_contrasts_over_all = mean(out.contrasts, 2);

figure;
plot(out.est.dithers, mean_contrasts_over_all, '-o');
title('Mean of the Mean Contrasts vs. Dithers in Voltages');
xlabel('Dither Voltage');
xlim([min(out.est.dithers)/2, max(out.est.dithers)*1.5]);
ylabel('Mean of Mean Contrasts');
grid on;

fprintf('Plots generated successfully.\n');

%% Save output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the current date and time to create a unique folder name
current_time = datestr(now, 'yyyy-mm-dd_HH-MM');
folder_name = ['dither_test_' current_time];
output_dir = fullfile(mp.path.config,'/', folder_name);

% Create the directory if it does not exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Save the generated plots
saveas(figure(1), fullfile(output_dir, 'Contrast_per_iteration.png'));
saveas(figure(2), fullfile(output_dir, 'DM1_Dither_std_vs_iteration.png'));
if any(mp.dm_ind == 2)
    saveas(figure(3), fullfile(output_dir, 'DM2_Dither_std_vs_iteration.png'));
    saveas(figure(4), fullfile(output_dir, 'Mean_Contrast_vs_Mean_DM2_Dither_std.png'));
    saveas(figure(5), fullfile(output_dir, 'Mean_of_Mean_Contrasts_vs_Dithers.png'));
else
    saveas(figure(3), fullfile(output_dir, 'Mean_Contrast_vs_Mean_DM1_Dither_std.png'));
    saveas(figure(4), fullfile(output_dir, 'Mean_of_Mean_Contrasts_vs_Dithers.png'));
end

% Save the MATLAB struct
save(fullfile(output_dir, 'mp_struct.mat'), 'mp');

fprintf('Plots, data, and MATLAB struct saved successfully in %s.\n', output_dir);

function [mp,ev] = pinned_act_safety_check(mp,ev)
ev.exit_flag = false;
% Update new pinned actuators
if any(mp.dm_ind == 1) 
    ev.dm1.new_pinned_actuators = setdiff(mp.dm1.pinned, ev.dm1.initial_pinned_actuators);
    ev.dm1.act_ele_pinned = mp.dm1.pinned(ismember(ev.dm1.new_pinned_actuators,mp.dm1.act_ele));
end
if any(mp.dm_ind == 2) 
    ev.dm2.new_pinned_actuators = setdiff(mp.dm2.pinned, ev.dm2.initial_pinned_actuators);
    ev.dm2.act_ele_pinned = mp.dm2.pinned(ismember(ev.dm2.new_pinned_actuators,mp.dm2.act_ele));

end


%  Check that no new actuators have been pinned
if size(ev.dm1.new_pinned_actuators,2)>0 || size(ev.dm2.new_pinned_actuators,2)>0

    % Print error warning
    if size(ev.dm1.new_pinned_actuators,2)>0
        fprintf('New DM1 pinned: [%s]\n', join(string(ev.dm1.new_pinned_actuators), ','));
    end
    
    if size(ev.dm2.new_pinned_actuators,2)>0
        fprintf('New DM2 pinned: [%s]\n', join(string(ev.dm2.new_pinned_actuators), ','));
    end

    % If new pinned actuators are used in jacobian, set flag to true so no image is taken 
    if size(ev.dm1.act_ele_pinned,2)>0 || size(ev.dm2.act_ele_pinned,2)>0
        save(fullfile([mp.path.config,'/','/ev_exit_',num2str(mp.dither),'.mat']),'ev')
        save(fullfile([mp.path.config,'/','/mp_exit_',num2str(mp.dither),'.mat']),"mp")

        % Reset pinned actuators to original
        if any(mp.dm_ind == 1)
            mp.dm1.pinned = ev.dm1.initial_pinned_actuators;
        end
        if any(mp.dm_ind == 2)
            mp.dm2.pinned = ev.dm2.initial_pinned_actuators;
        end

        % Set exit flag to true
        ev.exit_flag = true;
    end
end

end
