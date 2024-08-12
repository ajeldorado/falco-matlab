%% Replace path with the latest experiment .all file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("C:\Users\chris\OneDrive\Documents\FALCO_results\EKF_Series1_Trial1\EKF_Series1_Trial1_config.mat")
% load other files too
fprintf('Successfully loaded files.\n');


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
Itr = 3; 

%% Initializing necessary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subbandImage = falco_get_sbp_image(mp, iSubband)
    if mp.flagSim
        subbandImage = falco_get_sim_sbp_image(mp, iSubband);
    else
        subbandImage = falco_get_testbed_sbp_image(mp, iSubband);
    end
end

function heightMap = falco_calc_act_height_from_voltage(dm)
    switch lower(dm.fitType)
        case {'linear', 'poly1'}
            heightMap = dm.VtoH .* dm.V;
        case {'quadratic', 'poly2'}
            if ~isfield(dm, 'p1') || ~isfield(dm, 'p2') || ~isfield(dm, 'p3')
                error("The fields p1, p2, and p3 must exist when dm.fitType == 'quadratic'.\n" + ...
                    "Those fields satisfy the formula:\n" + ...
                    "height = p1*V*V + p2*V + p3");
            end
            Vbias = dm.biasMap;
            Vtotal = dm.V + Vbias;
            heightMapTotal = dm.p1.*Vtotal.^2 + dm.p2.*Vtotal + dm.p3;
            heightMapBias = dm.p1.*Vbias.^2 + dm.p2.*Vbias + dm.p3;
            heightMap = heightMapTotal - heightMapBias;
        case {'fourier2'}
            if ~isfield(dm, 'a0') || ~isfield(dm, 'a1') || ~isfield(dm, 'a2') || ...
                    ~isfield(dm, 'b1') || ~isfield(dm, 'b2') || ~isfield(dm, 'w')
                error("The fields a0, a1, a2, b1, b2, and w must exist when dm.fitType == 'fourier2'.\n" + ...
                    "Those fields satisfy the formula:\n" + ...
                    "height = a0 + a1*cos(V*w) + b1*sin(V*w) + a2*cos(2*V*w) + b2*sin(2*V*w)");
            end  
            Vbias = dm.biasMap;
            Vtotal = dm.V + Vbias;
            heightMapTotal = dm.a0 + dm.a1.*cos(Vtotal.*dm.w) + dm.b1.*sin(Vtotal.*dm.w) + dm.a2.*cos(2*Vtotal.*dm.w) + dm.b2.*sin(2*Vtotal.*dm.w);
            heightMapBias = dm.a0 + dm.a1.*cos(Vbias.*dm.w) + dm.b1.*sin(Vbias.*dm.w) + dm.a2.*cos(2*Vbias.*dm.w) + dm.b2.*sin(2*Vbias.*dm.w);
            heightMap = heightMapTotal - heightMapBias;        
        otherwise
            error('Value of dm.fitType not recognized.');
    end
end

fprintf('Necessary functions initialized.\n');

%% Setting initial DM commands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Access the 2D voltage maps for specified iteration
init_command_dm1 = out.dm1.Vall(:,:,Itr);  % 2D voltage map for dm1 at iteration Itr
init_command_dm2 = out.dm2.Vall(:,:,Itr);  % 2D voltage map for dm2 at iteration Itr

% Load the second to last command
delta_dm1 = abs(init_command_dm1-out.dm1.Vall(:,:,Itr-1));
delta_dm2 = abs(init_command_dm2-out.dm2.Vall(:,:,Itr-1));

mp.dm1.V_dz = init_command_dm1;
mp.dm2.V_dz = init_command_dm2;

% delta_dm_mean = mean([std(delta_dm1(init_command_dm1 ~= 0)), std(delta_dm2(init_command_dm2 ~= 0))]);
delta_dm_mean = mean([mean(delta_dm1(init_command_dm1 ~= 0)), mean(delta_dm2(init_command_dm2 ~= 0))]);

fprintf('Successfully loaded initial DM commands.\n');

%% Setting initial variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mp.est.dithers = delta_dm_mean*[1, 2, 3, 5, 10]; % ideally 3 dithers
% mp.est.dithers = [1e-10, 1e-5, 1e-3, 1];
num_dithers = numel(mp.est.dithers); % set to the same index number of mp.est.dithers
num_iterations = 100;
% num_subbands = 1;
iSubband = ceil(mp.Nsbp / 2); % Pick the middle wavelength
fprintf('Now setting necessary initial variables.\n');

%% Set total command for estimator image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: need to save these commands for each iteration separately

if Itr > 1
    if ~isfield(mp.dm1,'dV'); mp.dm1.dV = zeros(mp.dm1.Nact); end
    if ~isfield(mp.dm2,'dV'); mp.dm2.dV = zeros(mp.dm2.Nact); end
    % efc_command = get_dm_command_vector(mp,mp.dm1.dV, mp.dm2.dV);
else
    efc_command = 0 * mp.est.dithers(1);
    mp.dm1.dV = zeros(size(init_command_dm1));
    mp.dm2.dV = zeros(size(init_command_dm2));
end

fprintf('Now entering first for loop. (give statistics of each dither)\n');


%% Create an array for mp.output arrays
% contrasts --> # dithers (rows) vs # iterations (columns)

% this needs to be 2 dimensional
contrasts = zeros(num_dithers, num_iterations);
dithers_nm_1 = zeros(num_dithers, num_iterations);
dithers_nm_2 = zeros(num_dithers, num_iterations);

%% For loops %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One loop that goes over the iterations
for i = 1:num_dithers
    dither = mp.est.dithers(i); % assuming only one dither value, otherwise use i if there are multiple dithers

    % Second loop that will give statistics of each dither
    for j = 1:num_iterations
        fprintf('Now inside second for loop.\n');
        
        V1 = zeros(mp.dm1.NactTotal,1);
        V2 = zeros(mp.dm2.NactTotal,1);

        % Calculate the dither commands
        if any(mp.dm_ind == 1)
            V1(mp.dm1.act_ele) = normrnd(0,mp.est.dithers(i),size(mp.dm1.act_ele));
            DM1Vdither = reshape(V1,[mp.dm1.Nact,mp.dm1.Nact]);
        else
            DM1Vdither = zeros(size(mp.dm1.V));
        end

        if any(mp.dm_ind == 2)
            V2(mp.dm2.act_ele) = normrnd(0,mp.est.dithers(i),size(mp.dm2.act_ele));
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
            I = falco_get_sbp_image(mp, iSubband); % iSubband : index of subband for which to take the image
            disp('falco_get_sbp_image called successfully');

            % Need to convert if it's not in contrast units
            % I0 = ev.imageArray / peak to put it into contrast units % ev is the estimator variable... (internal)
            mean_contrast = mean(I(mp.Fend.score.mask)); % this is the dark zone mask
            disp(['Mean contrast for subband ', num2str(iSubband), ': ', num2str(mean_contrast)]);
            
            % Generate command to apply to DMs
            if any(mp.dm_ind == 1)
                % Note falco_set_constrained_voltage does not apply the command to the DM
                mp.dm1 = falco_set_constrained_voltage(mp.dm1, DM1Vdither); 
            end
    
            if any(mp.dm_ind == 2)
                mp.dm2 = falco_set_constrained_voltage(mp.dm2, DM2Vdither); 
            end
            
            %% STILL NEED TO DO:
            DM1Vdither_nm = falco_calc_act_height_from_voltage(mp.dm1)*10^9; % convert to nano meters
            DM2Vdither_nm = falco_calc_act_height_from_voltage(mp.dm2)*10^9;

            % Store results
            contrasts(i, j) = mean_contrast;

            % Calculate standard deviation of the dithered commands (excluding zeros)
            dithers_nm_1(i, j) = std(DM1Vdither_nm(init_command_dm1 ~= 0));
            dithers_nm_2(i, j) = std(DM2Vdither_nm(init_command_dm2 ~= 0));
        end
        
    end
end

%% Plotting
% four things:
%   1. contrast per iteration for each of the dithers
figure;
for i = 1:num_dithers
    semilogy(1:num_iterations, contrasts(i, :), '-o');
    hold on;
end
title('Contrast per iteration for each of the dithers');
xlabel('Iteration');
ylabel('Contrast');
legend(arrayfun(@(x) ['Dither ' num2str(x)], 1:num_dithers, 'UniformOutput', false));
hold off;

%   2. dither std nms vs. iteration for each of the dithers
figure;
for i = 1:num_dithers
    plot(1:num_iterations, dithers_nm_1(i, :), '-o');
    hold on;
end
title('DM1 Dither std nms vs. iteration for each of the dithers');
xlabel('Iteration');
ylabel('Standard Deviation (nm)');
legend(arrayfun(@(x) ['Dither ' num2str(x)], 1:num_dithers, 'UniformOutput', false));
hold off;

figure;
for i = 1:num_dithers
    plot(1:num_iterations, dithers_nm_2(i, :), '-o');
    hold on;
end
title('DM2 Dither std nms vs. iteration for each of the dithers');
xlabel('Iteration');
ylabel('Standard Deviation (nm)');
legend(arrayfun(@(x) ['Dither ' num2str(x)], 1:num_dithers, 'UniformOutput', false));
hold off;

%   3. take the averages along the dimensions of both the contrasts and the dither std
%                   nms over the iteration axis, then you plot the mean std
%                   dithers (as x-axis) and mean contrast (as y-axis)
mean_contrasts = mean(contrasts, 2);
mean_dithers_nm_1 = mean(dithers_nm_1, 2);
mean_dithers_nm_2 = mean(dithers_nm_2, 2);

figure;
semilogy(mean_dithers_nm_1, mean_contrasts, '-o');
title('Mean Contrast vs. Mean DM1 Dither std nms');
xlabel('Mean DM1 Dither std (nm)');
xlim([min(mean_dithers_nm_1)/2, max(mean_dithers_nm_1)*1.5]);
ylabel('Mean Contrast');
grid on;

figure;
plot(mean_dithers_nm_2, mean_contrasts, '-o');
title('Mean Contrast vs. Mean DM2 Dither std nms');
xlabel('Mean DM2 Dither std (nm)');
xlim([min(mean_dithers_nm_2)/2, max(mean_dithers_nm_2)*1.5]);
ylabel('Mean Contrast');
grid on;

%   4. mean of the mean contrasts (y-axis) vs. dithers in voltages (x-axis)
mean_contrasts_over_all = mean(contrasts, 2);

figure;
semilogy(mp.est.dithers, mean_contrasts_over_all, '-o');
title('Mean of the Mean Contrasts vs. Dithers in Voltages');
xlabel('Dither Voltage');
xlim([min(mp.est.dithers)/2, max(mp.est.dithers)*1.5]);
ylabel('Mean of Mean Contrasts');
grid on;

fprintf('Plots generated successfully.\n');



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