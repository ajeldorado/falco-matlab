% Things to check:
% EFC command
% Pinned actuators
% Drift command
% OL images (delta command between OL and CL)
% check if OL or CL estimate gets handed to controller

%%

data_path = ['C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55'];
load('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\EKF_Series8_Trial55_config.mat')
load('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\EKF_Series8_Trial55_snippet.mat')


%%
num_iter = 420;
iter_arr = [1:num_iter].';
weird_iter = [317,375,403,417]-1;

i_bb = out.IrawScoreHist(1:num_iter);

i_bb_est = mean(out.normIntMeasScore,2); % double check if this is OL or CL estimate
i_bb_est = i_bb_est(1:num_iter);

figure(1)
semilogy(iter_arr,i_bb,iter_arr,i_bb_est)
legend('measured','estimated','location','northwest')
title('Measured and Estimated Mean DZ contrast')
xlabel('iteration')
ylabel('contrast')

figure(2)

semilogy(iter_arr(weird_iter(1)-2:weird_iter(1)+2),i_bb(weird_iter(1)-2:weird_iter(1)+2), ...
    iter_arr(weird_iter(1)-2:weird_iter(1)+2),i_bb_est(weird_iter(1)-2:weird_iter(1)+2))
legend('measured','estimated','location','northwest')
title('Measured and Estimated Mean DZ contrast Zoom')
xlabel('iteration')
ylabel('contrast')

% Based on these plots it does not seem to be the estimate that is driving
% the divergence



%%
probe_iter = ['probing_data_',num2str(weird_iter(end)+1),'.mat'];
probe_iter_path = fullfile('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\',probe_iter);
load(probe_iter_path)

i_bb_ol = mean(ev.IOLScoreHist,2);
i_bb_ol = i_bb_ol(1:weird_iter(end));


figure(10)
semilogy(iter_arr,i_bb,iter_arr(1:weird_iter(end)),i_bb_ol,iter_arr,i_bb_est)
legend('measured CL','measured OL','estimated','location','northwest')
title('Measured and Estimated Mean DZ contrast')
xlabel('iteration')
ylabel('contrast')
% xlim([weird_iter(1)-2,weird_iter(1)+2])

% Based on this figure it does not *appear* to be related to taking OL
% images

%% pinned act check


num_pinned_act = out.dm1.Npinned(1:num_iter);

figure(100)
yyaxis left
semilogy(iter_arr,i_bb)
ylabel('contrast')

yyaxis right
plot(iter_arr,num_pinned_act)
ylabel('number of actuators pinned')
ylim([min(num_pinned_act),max(num_pinned_act)+30])
legend('contrast','pinned actuators')
title('Mean DZ contrast and pinned actuators')
xlabel('iteration')

%%
wi = 4;
pinned_actm1 = cell2mat(out.dm1.pinned(weird_iter(wi)-1));
pinned_act1 = cell2mat(out.dm1.pinned(weird_iter(wi)));
pinned_actp1 = cell2mat(out.dm1.pinned(weird_iter(wi)+1));
dm_arr_zero = zeros(size(mp.dm1.V));

dm_arr1 = dm_arr_zero;
dm_arr1(pinned_act1) = 1;

dm_arrp1 = dm_arr_zero;
dm_arrp1(pinned_actp1) = 1;

dm_arrm1 = dm_arr_zero;
dm_arrm1(pinned_actm1) = 1;

figure(111)
subplot(131)
imagesc(dm_arrm1)
title(num2str(weird_iter(wi)-1))
subplot(132)
imagesc(dm_arr1)
title([{'pinned act'},{num2str(weird_iter(wi))}])
subplot(133)
imagesc(dm_arrp1)
title(num2str(weird_iter(wi)+1))

figure(112)
subplot(121)
imagesc(dm_arrm1-dm_arr1)
title([{'new pinned act'},{[num2str(weird_iter(wi)-1),'-',num2str(weird_iter(wi))]}])
subplot(122)
imagesc(dm_arr1-dm_arrp1)
title([{'new pinned act'},{[num2str(weird_iter(wi)),'-',num2str(weird_iter(wi)+1)]}])


% Pontentially a correlation between pinned actuators and steps


%% New pinned
start_check = 300;
dm_new = dm_arr_zero;
figure(222)
imagesc(dm_new)
ind = 1;
pinned0 = cell2mat(out.dm1.pinned(1));
new_pinned = zeros(30,num_iter-start_check);

new_pinned_iter = zeros(30,num_iter-start_check);
new_pinned_volts_iter = zeros(30,num_iter-start_check);
temp0 = pinned0;
for i = start_check:1:num_iter
    temp = cell2mat(out.dm1.pinned(i));
    temp_volt = cell2mat(out.dm1.Vpinned(i));

    new_temp = setdiff(temp,pinned0);
    new_pinned(1:size(new_temp,1),ind) = new_temp;

    [new_temp_iter,ind_temp] = setdiff(temp,temp0);
    
    new_pinned_iter(1:size(new_temp_iter,1),ind) = setdiff(temp,temp0);
    new_pinned_volts_iter(1:size(new_temp_iter,1),ind) = temp_volt(ind_temp); 

    dm_new(new_temp) = 1;
    imagesc(dm_new)
    title(num2str(i))
%     drawnow
    pause(0.1)
    ind = ind+1;
    temp0 = temp;
end

%% load mean contrast from divergence_pinned_act_check

for i = 1:1:num_iter-start_check

temp_acts = new_pinned_iter(:,i);temp_acts = temp_acts(temp_acts~=0);

temp_volts = new_pinned_volts_iter(:,i);temp_volts = temp_volts(temp_volts~=0);

% try
    temp_acts_pmax = temp_acts(temp_volts>0);
    temp_acts_pmin = temp_acts(temp_volts<0);
    if size(temp_acts_pmax,1)>0
        [a,inds_max] = union(pinned_act_end , temp_acts_pmax);
        contrast_pmax = mean_contrast_pmax(inds_max);
    else
        contrast_pmax = 0;
    end
    if size(temp_acts_pmin,1)>0
        [a,inds_min] = union(pinned_act_end , temp_acts_pmin);
        contrast_pmin = mean_contrast(inds_min);
    else
        contrast_pmin = 1e-8;
    end

    contrast_iter(i,1) = sum([contrast_pmax;contrast_pmin]);
% catch
%     contrast_iter(i,1) = 1e-8;
% end


end
%%
red_array = start_check-1:num_iter-2;
contrast_iter_steps_only = contrast_iter;
contrast_iter_steps_only(contrast_iter_steps_only<1.1e-8) = 0;
contrast_pred = zeros(num_iter,1);
contrast_pred(1:start_check,1) = i_bb(1:start_check,1);  
contrast_pred(red_array,1) = i_bb(red_array,1)+contrast_iter_steps_only;  
figure;
semilogy(iter_arr,contrast_pred,'r')
hold on
semilogy(iter_arr,i_bb,'k')
semilogy(red_array(contrast_iter>1.1e-8),contrast_iter(contrast_iter>1.1e-8),'r*')
hold off
xlabel('iteration')
ylabel('mean dz contrast')
legend('predicted from sims based on pinned act','experiment')
title('Predicted change in contrast based on pinned actuators')
%%

figure
hold on
for i = 1:1:num_iter-start_check
    temp = new_pinned_iter(:,i);

    try
    plot(i+start_check,temp(temp~=0),'k*')
    catch
    end

end
hold off

figure
hold on
for i = 1:1:num_iter-start_check
    temp = new_pinned_volts_iter(:,i);

    try
    plot(i+start_check,temp(temp~=0),'k*')
    catch
    end

end
hold off


%% drift check
wi = 4;

drift_witer1_name = ['drift_command_it',num2str(weird_iter(wi)),'.fits'];
drift_witer1_path = fullfile('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\',drift_witer1_name);

drift_witer1m1_name = ['drift_command_it',num2str(weird_iter(wi)-1),'.fits'];
drift_witer1m1_path = fullfile('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\',drift_witer1m1_name);


drift_witer1 = fitsread(drift_witer1_path);
drift_witer1m1 = fitsread(drift_witer1m1_path);

drift_witer1 = drift_witer1(:,:,1);
drift_witer1m1 = drift_witer1m1(:,:,1);

figure(101)
plot(drift_witer1(:),'*')
hold on
plot(drift_witer1m1(:))
plot((drift_witer1m1(:)) - drift_witer1(:))
hold off
xlabel('actuator')
ylabel('voltage')
legend(num2str(weird_iter(wi)),num2str(weird_iter(wi)-1))


% sanity check that drift commands are different
if wi == 4
    drift_witer1p1_name = ['drift_command_it',num2str(weird_iter(wi)+1),'.fits'];
    drift_witer1p1_path = fullfile('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\',drift_witer1p1_name);
    drift_witer1p1 = fitsread(drift_witer1p1_path);
    drift_witer1p1 = drift_witer1p1(:,:,1);

    d_drift_m11 = drift_witer1m1(:) - drift_witer1(:);
    d_drift_1p1 = drift_witer1(:) - drift_witer1p1(:);

    figure(102)
    plot(d_drift_m11,'*')
    hold on
    plot(d_drift_1p1)
    hold off
    xlabel('actuator')
    ylabel('voltage')
    legend('m11','1p1')
end


% Drift command seems ok

%% EFC command check

wi = 4;
wi_pm1 = weird_iter(wi)-1;

efc_witer1_name = ['efc_command_it',num2str(weird_iter(wi)),'.fits'];
efc_witer1_path = fullfile('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\',efc_witer1_name);

efc_witer1m1_name = ['efc_command_it',num2str(wi_pm1),'.fits'];
efc_witer1m1_path = fullfile('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\',efc_witer1m1_name);


efc_witer1 = fitsread(efc_witer1_path);
efc_witer1m1 = fitsread(efc_witer1m1_path);

efc_witer1 = efc_witer1(:,:,1);
efc_witer1m1 = efc_witer1m1(:,:,1);

figure(101)
plot(efc_witer1(:),'*')
hold on
plot(efc_witer1m1(:))
plot((efc_witer1m1(:)) - efc_witer1(:))
hold off
xlabel('actuator')
ylabel('voltage')
legend(num2str(weird_iter(wi)),num2str(wi_pm1))

% efc command seems ok

%% dm1_V = V0 + dither + drift check

wi = 4;
wi_pm1 = weird_iter(wi)-1;

dm1_V__witer1_name = ['dm1_V_it',num2str(weird_iter(wi)),'.fits'];
dm1_V__witer1_path = fullfile('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\',dm1_V__witer1_name);

dm1_V__witer1m1_name = ['dm1_V_it',num2str(wi_pm1),'.fits'];
dm1_V__witer1m1_path = fullfile('C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55\',dm1_V__witer1m1_name);


dm1_V__witer1 = fitsread(dm1_V__witer1_path);
dm1_V__witer1m1 = fitsread(dm1_V__witer1m1_path);

dm1_V__witer1 = dm1_V__witer1(:,:,1);
dm1_V__witer1m1 = dm1_V__witer1m1(:,:,1);

figure(103)
plot(dm1_V__witer1(:),'*')
hold on
plot(dm1_V__witer1m1(:))
plot((dm1_V__witer1m1(:)) - dm1_V__witer1(:))
hold off
xlabel('actuator')
ylabel('voltage')
legend(num2str(weird_iter(wi)),num2str(wi_pm1),'diff')

% dm1_V looks fine

%% pinned highlighted full DM command
v_tot_m1 = efc_witer1m1 + dm1_V__witer1m1;
v_tot = efc_witer1 + dm1_V__witer1;

pinned_act_end = new_pinned(:,end);
pinned_act_end = pinned_act_end(pinned_act_end~=0);
figure(222)
subplot(121)
plot(v_tot_m1(pinned_act_end))
title('New pinned actuators during DZM iter 415')
xticks(1:length(pinned_act_end));
xticklabels(string(pinned_act_end))
xlabel('pinned actuator')
ylabel('voltage')
subplot(122)
plot(v_tot(pinned_act_end))
title('New pinned actuators during DZM iter 416')
xticks(1:length(pinned_act_end));
xticklabels(string(pinned_act_end))
xlabel('pinned actuator')
ylabel('voltage')

%%
% v_tot = efc_witer1 + dm1_V__witer1;
v_tot(pinned_act1) = NaN;

figure(104)
% imagesc(v_tot)
imagesc(v_tot,'AlphaData',~isnan(v_tot))
c = colorbar;
c.Label.String = 'voltage';
clim([-0.15,0.15])
title('Full DM command with pinned actuators in white')

%%

v_tot = efc_witer1 + dm1_V__witer1;
v_tot(pinned_act1) = NaN;

figure(104)
% imagesc(v_tot)
imagesc(v_tot,'AlphaData',~isnan(v_tot))
c = colorbar;
c.Label.String = 'voltage';
clim([-0.15,0.15])
title('Full DM command with pinned actuators in white')






%% Temp / humidity?
