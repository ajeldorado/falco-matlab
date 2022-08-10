% run a couple iterations of run_falco_efc_vortex_v1 first

load('new_pinned_trial55.mat')

pinned_act_end = new_pinned(:,end);
pinned_act_end = pinned_act_end(pinned_act_end~=0);

num_check = size(pinned_act_end,1);

initial_command = out.dm1.Vall(:,:,end);

imageArray_before = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1, mp.Nsbp);
mean_contrast_before = zeros(1, mp.Nsbp);

for iSubband = 1:mp.Nsbp
    imageArray_before(:,:,1,iSubband) = falco_get_sbp_image(mp, iSubband);
    I0 = imageArray_before(:,:,1,iSubband);
    mean_contrast_before(1,iSubband) = mean(I0(mp.Fend.score.mask));
end


imageArray = zeros(mp.Fend.Neta, mp.Fend.Nxi, num_check, mp.Nsbp);

mean_contrast = zeros(num_check, mp.Nsbp);

for i = 1:1:num_check
    mp.dm1.V = initial_command;
    mp.dm1.V(pinned_act_end(i)) = -0.02;
    
for iSubband = 1:mp.Nsbp
    imageArray(:,:,i,iSubband) = falco_get_sbp_image(mp, iSubband);
    I0 = imageArray(:,:,i,iSubband);
    mean_contrast(i,iSubband) = mean(I0(mp.Fend.score.mask));
end

end


imageArray_pmax = zeros(mp.Fend.Neta, mp.Fend.Nxi, num_check, mp.Nsbp);

mean_contrast_pmax = zeros(num_check, mp.Nsbp);

for i = 1:1:num_check
    mp.dm1.V = initial_command;
    mp.dm1.V(pinned_act_end(i)) = 0.1;
    
for iSubband = 1:mp.Nsbp
    imageArray_pmax(:,:,i,iSubband) = falco_get_sbp_image(mp, iSubband);
    I0 = imageArray_pmax(:,:,i,iSubband);
    mean_contrast_pmax(i,iSubband) = mean(I0(mp.Fend.score.mask));
end

end

%% Note which actuators have largest influence on contrast
sens_max_acts = pinned_act_end(mean_contrast_pmax>1e-8);
sens_min_acts = pinned_act_end(mean_contrast>1e-8);


fprintf('sensitive high actuators: [%s]\n', join(string(sens_max_acts), ','));
fprintf('sensitive low actuators: [%s]\n', join(string(sens_min_acts), ','));
fprintf('sensitive both actuators: [%s]\n', join(string(intersect(sens_min_acts,sens_max_acts)), ','));


%%
figure
subplot(121)
semilogy(1:num_check,mean_contrast_before(1,1)*ones(num_check,1),'c','linewidth',3)
hold on
semilogy(1:num_check,mean_contrast(:,1),'k*')
hold off

xticks(1:length(pinned_act_end));
xticklabels(string(pinned_act_end))
xlabel('pinned actuator')
ylabel('mean contrast')
title([{'mean contrast comparison for various'},{'pinned-to-ZERO actuators (SIM)'}])
legend('initial contrast','contrast with one pinned actuator')

subplot(122)

semilogy(1:num_check,mean_contrast(:,1) - mean_contrast_before(1,1),'*')
xticks(1:length(pinned_act_end));
xticklabels(string(pinned_act_end))
xlabel('pinned actuator')
ylabel('change in mean contrast (SIM)')
title([{'c_{pinned}-c_0'},{'pinned-to-ZERO actuators (SIM)'}])


print()
%%
figure
subplot(121)
semilogy(1:num_check,mean_contrast_before(1,1)*ones(num_check,1),'c','linewidth',3)
hold on
semilogy(1:num_check,mean_contrast_pmax(:,1),'k*')
hold off

xticks(1:length(pinned_act_end));
xticklabels(string(pinned_act_end))
xlabel('pinned actuator')
ylabel('mean contrast')
title([{'mean contrast comparison for various'},{'pinned-to-MAX actuators (SIM)'}])
legend('initial contrast','contrast with one pinned actuator')

subplot(122)

semilogy(1:num_check,mean_contrast_pmax(:,1) - mean_contrast_before(1,1),'*')
xticks(1:length(pinned_act_end));
xticklabels(string(pinned_act_end))
xlabel('pinned actuator')
ylabel('change in mean contrast')
title([{'c_{pinned}-c_0'},{'pinned-to-MAX actuators (SIM)'}])

%%
dm_temp = 0*initial_command;
dm_temp(sens_max_acts) = 1;
dm_temp(sens_min_acts) = -1;
figure;
imagesc(dm_temp)
title([{'Acts with large influence when pinned'},{'1=pinned high, -1=pinned high or low'}])
colorbar

