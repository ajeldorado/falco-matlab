%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test dithers to find contrast-dither relationship on IACT.
% Generate dark zone using preferred method and then run this script.
%
% sfr 01082022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up data paths
dhmOutDir = '/proj/iact/data/dhm/driftTests/';
runTag = datestr(datetime('now'),30); % get a string with the current date and time 


%% Get mean DM command between last two iterations
dV1_mean = mean(mean(abs(out.dm1.Vall(:,:,end) - out.dm1.Vall(:,:,end-1))));
dV2_mean = mean(mean(abs(out.dm2.Vall(:,:,end) - out.dm2.Vall(:,:,end-1))));
dV_mean = mean([dV1_mean,dV2_mean]);

% p1 = 5.587e-6;
% p2 = 0.14088e-6;
% dx_mean = 2*p1*dV1_mean + p2;
% disp("mean displacement of DM actuators between last two iterations",num2str(dx_mean))


%%
% dither_1 = [0, 0.5*dV1_mean:dV1_mean:dV1_mean*10];
% dither_2 = [0, 0.5*dV2_mean:dV2_mean:dV2_mean*10];

dither = [0, 0.5*dV_mean:dV_mean:dV_mean*10];

initial_command_1 = out.dm1.Vall(:,:,end);
initial_command_2 = out.dm2.Vall(:,:,end);


% Pick the middle wavelength
si = ceil(mp.Nsbp / 2);

active_dm_arr = [1,0; 0,1; 1,1]; % needs to be same size as number of DMs, 1 means DMi is dithering
contrast = zeros(50,size(dither,1),3);
mean_contrasts = zeros(1,size(dither,1),3);
for i = 1:1:3
    active_dm = active_dm_arr(i,:);
    %     if sum(active_dm) == 2
    %         dither = mean([dither_1;dither_2]);
    %         dither_1 = dither;
    %         dither_2 = dither;
    %     end
    figure(i)
    contrast(:,:,i) = get_contrast_vs_dither(mp,dither,initial_command_1,initial_command_2,si,active_dm);

    mean_contrasts(:,:,i) = mean(contrast(:,2:end,i),1);


end
%% Plot final result for DM1

dither_legend = "dither = " + string(dither);

figure(2)
subplot(1,2,1)
semilogy([1:1:50],contrast(:,:,1))
xlim([1,50])
legend(dither_legend)
title([{'DM1'},{'Mean dz contrast vs image number'},{'for varying dither'}])
xlabel('image number')
ylabel('mean dz contrast')


subplot(1,2,2)
plot(dither,mean_contrasts_1(:,:,1))
xlabel('$\sigma_{dither} [V/\sqrt{iter}]$','Interpreter','latex')
ylabel('mean dz contrast')
title([{'DM1'},{'mean dz contrast'},{'vs dither std dev'}])

%% Save
out_dither_test_1.contrast = contrast(:,:,1);
out_dither_test_1.mean_contrast = mean_contrasts(:,:,1);
out_dither_test_1.dither = dither;

out_dither_test_2.contrast = contrast(:,:,2);
out_dither_test_2.mean_contrast = mean_contrasts(:,:,2);
out_dither_test_2.dither = dither;

out_dither_test_12.contrast = contrast(:,:,3);
out_dither_test_12.mean_contrast = mean_contrasts(:,:,3);
out_dither_test_12.dither = dither;

save(fullfile(dhmOutDir,['out_dither_test_1',runTag,'.mat']),'out_dither_test_1')
save(fullfile(dhmOutDir,['out_dither_test_2',runTag,'.mat']),'out_dither_test_2')
save(fullfile(dhmOutDir,['out_dither_test_12',runTag,'.mat']),'out_dither_test_12')


%% Function to get contrast vs dither for a particular DM
function contrast = get_contrast_vs_dither(mp,dither,initial_command_1,initial_command_2,si,active_dm)

dm1_gain = active_dm(1);
dm2_gain = active_dm(2);

contrast = zeros(50,size(dither,1));
for i = 1:1:size(dither,1)
    
    for j = 1:1:50
        % Initialize voltage vector
        V1 = zeros(mp.dm1.NactTotal,1);
        V2 = zeros(mp.dm2.NactTotal,1);

        % set random command for active actuators?
        V1(mp.dm1.act_ele) = dm1_gain * normrnd(0,dither(i),size(mp.dm1.act_ele));
        V2(mp.dm2.act_ele) = dm2_gain * normrnd(0,dither(i),size(mp.dm2.act_ele));

        % Update mp.dm1 object
        mp.dm1.V = initial_command_1 + reshape(V1,[mp.dm1.Nact,mp.dm1.Nact]);
        mp.dm2.V = initial_command_2 + reshape(V2,[mp.dm2.Nact,mp.dm2.Nact]);

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
        xlim([1,50])
        ylim([0.1*min(contrast(:,i)),2*max(contrast(:,i))])
        xlabel('image number')
        ylabel('mean dark zone contrast')
        title([{'mean dark zone contrast'}, {'for $\sigma_{dither}$ = ',num2str(dither(i))}],'Interpreter','latex')

        drawnow
    end


end
end
