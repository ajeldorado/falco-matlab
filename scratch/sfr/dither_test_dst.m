%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test dithers to find contrast-dither relationship on IACT.
% Generate dark zone using preferred method and then run this script.
%
% sfr 01082022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up data paths

dhmOutDir = fullfile(mp.path.config,'dither_test',mp.runLabel);
runTag = datestr(datetime('now'),30); % get a string with the current date and time 

mkdir(dhmOutDir)

%% Get active actuators for DM1 and DM2
% mp.dm_ind = [1,2];
% % Normalization and throughput calculations
% 
% mp = falco_compute_psf_norm_factor(mp);
% 
% [mp, thput, ImSimOffaxis] = falco_compute_thput(mp);
% % mp.thput_vec(Itr) = max(thput); % note: max() needed when mp.flagFiber==true
% 
% % Control Jacobian
% 
% mp = falco_set_jacobian_modal_weights(mp);
% jacStruct =  model_Jacobian(mp);
% 
% 
% [mp, jacStruct] = falco_ctrl_cull_weak_actuators(mp, cvar, jacStruct);
%% Initialize pinned act list

% Initialize pinned actuator check
ev.dm1.initial_pinned_actuators = mp.dm1.pinned;
if any(mp.dm_ind == 2); ev.dm2.initial_pinned_actuators = mp.dm2.pinned; end

ev.dm1.new_pinned_actuators = [];
ev.dm2.new_pinned_actuators = [];
ev.dm1.act_ele_pinned = [];
ev.dm2.act_ele_pinned = [];

%% Get mean DM command between last two iterations
dV1_mean = mean(mean(abs(out.dm1.Vall(:,:,end) - out.dm1.Vall(:,:,end-1))));
if any(mp.dm_ind == 2)
    dV2_mean = mean(mean(abs(out.dm2.Vall(:,:,end) - out.dm2.Vall(:,:,end-1)))); 
else
    dV2_mean = dV1_mean;
end
dV_mean = mean([dV1_mean,dV2_mean]);

% p1 = 5.587e-6;
% p2 = 0.14088e-6;
% dx_mean = 2*p1*dV1_mean + p2;
% disp("mean displacement of DM actuators between last two iterations",num2str(dx_mean))


%%

% Set number of exposures to average over
num_iter = 5;

dither = [0, 0.5*dV_mean:dV_mean:dV_mean*10].';

initial_command_1 = out.dm1.Vall(:,:,end);
initial_command_2 = out.dm2.Vall(:,:,end);


% Pick the middle wavelength
si = ceil(mp.Nsbp / 2);

active_dm_arr = [1,0; 0,1; 1,1]; % needs to be same size as number of DMs, 1 means DMi is dithering
contrast = zeros(num_iter,size(dither,1),size(active_dm_arr,1));
mean_contrasts = zeros(1,size(dither,1),size(active_dm_arr,1));
for i = 1:1:size(active_dm_arr,1)
    active_dm = active_dm_arr(i,:);

    figure(i)
    contrast(:,:,i) = get_contrast_vs_dither(mp,dither,initial_command_1,initial_command_2,si,active_dm,num_iter,ev);

%     mean_contrasts(:,:,i) = mean(contrast(:,2:end,i),1);
    mean_contrasts(:,:,i) = mean(contrast(:,:,i),1);


end

%% Save
out_dither_test_1.contrast = contrast(:,:,1);
out_dither_test_1.mean_contrast = mean_contrasts(:,:,1);
out_dither_test_1.dither = dither;
save(fullfile(dhmOutDir,['out_dither_test_1',runTag,'.mat']),'out_dither_test_1')

out_dither_test_2.contrast = contrast(:,:,2);
out_dither_test_2.mean_contrast = mean_contrasts(:,:,2);
out_dither_test_2.dither = dither;
save(fullfile(dhmOutDir,['out_dither_test_2',runTag,'.mat']),'out_dither_test_2')

out_dither_test_12.contrast = contrast(:,:,3);
out_dither_test_12.mean_contrast = mean_contrasts(:,:,3);
out_dither_test_12.dither = dither;
save(fullfile(dhmOutDir,['out_dither_test_12',runTag,'.mat']),'out_dither_test_12')


%% Plot final result for DM1

dither_legend = "dither = " + string(dither);
figure(1)
subplot(1,2,1)
semilogy([1:1:num_iter],contrast(:,:,1))
xlim([1,num_iter])
legend(dither_legend)
title([{'DM1'},{'Mean dz contrast vs image number'},{'for varying dither'}])
xlabel('image number')
ylabel('mean dz contrast')


subplot(1,2,2)
plot(dither,mean_contrasts(:,:,1))
xlabel('$\sigma_{dither} [V/\sqrt{iter}]$','Interpreter','latex')
ylabel('mean dz contrast')
title([{'DM1'},{'mean dz contrast'},{'vs dither std dev'}])


figure(2)
subplot(1,2,1)
semilogy([1:1:num_iter],contrast(:,:,2))
xlim([1,num_iter])
legend(dither_legend)
title([{'DM2'},{'Mean dz contrast vs image number'},{'for varying dither'}])
xlabel('image number')
ylabel('mean dz contrast')


subplot(1,2,2)
plot(dither,mean_contrasts(:,:,2))
xlabel('$\sigma_{dither} [V/\sqrt{iter}]$','Interpreter','latex')
ylabel('mean dz contrast')
title([{'DM2'},{'mean dz contrast'},{'vs dither std dev'}])

figure(3)
subplot(1,2,1)
semilogy([1:1:num_iter],contrast(:,:,3))
xlim([1,num_iter])
legend(dither_legend)
title([{'DM12'},{'Mean dz contrast vs image number'},{'for varying dither'}])
xlabel('image number')
ylabel('mean dz contrast')


subplot(1,2,2)
plot(dither,mean_contrasts(:,:,3))
xlabel('$\sigma_{dither} [V/\sqrt{iter}]$','Interpreter','latex')
ylabel('mean dz contrast')
title([{'DM12'},{'mean dz contrast'},{'vs dither std dev'}])


%% Function to get contrast vs dither for a particular DM
function contrast = get_contrast_vs_dither(mp,dither,initial_command_1,initial_command_2,si,active_dm,num_iter,ev)

dm1_gain = active_dm(1);
dm2_gain = active_dm(2);

contrast = zeros(num_iter,size(dither,1));
for i = 1:1:size(dither,1)
    
    for j = 1:1:num_iter
        % Initialize voltage vector
        V1 = zeros(mp.dm1.NactTotal,1);
        V2 = zeros(mp.dm2.NactTotal,1);

        % set random command for active actuators?
        V1(mp.dm1.act_ele) = dm1_gain * normrnd(0,dither(i),size(mp.dm1.act_ele));
        V2(mp.dm2.act_ele) = dm2_gain * normrnd(0,dither(i),size(mp.dm2.act_ele));

        % Update mp.dm1 object
        mp.dm1.V = initial_command_1 + reshape(V1,[mp.dm1.Nact,mp.dm1.Nact]);
        mp.dm2.V = initial_command_2 + reshape(V2,[mp.dm2.Nact,mp.dm2.Nact]);
        
        % Check pinned actuators
        mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V); 
        mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V); 
        
        mp.dither = dither(i);
        [mp, ev] = pinned_act_safety_check(mp,ev);
        % Only take image if no actuators in act_ele will get pinned
        if ~ev.exit_flag
            % Take image with random command
            image = falco_get_sbp_image(mp,si); % this returns image in contrast units
            contrast(j,i) = mean(image(mp.Fend.score.mask));
            
            subplot(1,2,1)
            imshow(log10(abs(image)))
            title(['image for $\sigma_{dither}$ = ',num2str(dither(i))],'Interpreter','latex')
            clim([min(min(log10(abs(image)))), max(max(log10(abs(image))))])
            colorbar

            subplot(1,2,2)
            semilogy(1:j,contrast(1:j,i))
            xlim([1,num_iter])
            ylim([0.1*min(contrast(:,i)),2*max(contrast(:,i))])
            xlabel('image number')
            ylabel('mean dark zone contrast')
            title([{'mean dark zone contrast'}, {'for $\sigma_{dither}$ = ',num2str(dither(i))}],'Interpreter','latex')

            drawnow
        else    
            contrast(j,i) = 0;%contrast(j-1,i);%nan;
        end   
        


    end

    


end
end

function [mp,ev] = pinned_act_safety_check(mp,ev)
ev.exit_flag = false;
% Update new pinned actuators
if any(mp.dm_ind == 1) || any(mp.dm_drift_ind == 1)
    ev.dm1.new_pinned_actuators = setdiff(mp.dm1.pinned, ev.dm1.initial_pinned_actuators);
    ev.dm1.act_ele_pinned = mp.dm1.pinned(ismember(ev.dm1.new_pinned_actuators,mp.dm1.act_ele));
end
if any(mp.dm_ind == 2) || any(mp.dm_drift_ind == 2)
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
        mp.dm1.pinned = ev.dm1.initial_pinned_actuators;
        mp.dm2.pinned = ev.dm2.initial_pinned_actuators;
        
        % Set exit flag to true
        ev.exit_flag = true;
    end
end

end


