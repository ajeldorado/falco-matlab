%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test dithers to find contrast-dither relationship on IACT.
% Generate dark zone using preferred method and then run this script.
%
% sfr 01082022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dhmOutDir = '/proj/iact/data/dhm/driftTests/';
runTag = datestr(datetime('now'),30); % get a string with the current date and time 

dV_mean = mean(mean(abs(out.dm1.Vall(:,:,end) - out.dm1.Vall(:,:,end-1))));
p1 = 5.587e-6;
p2 = 0.14088e-6;
dx_mean = 2*p1*dV_mean + p2;
% GR suggests dx_mean = mp.dm1.VtoH*dV_mean;
disp("mean displacement of DM actuators between last two iterations",num2str(dx_mean))
%%
n_iter = 200;

figure(231)

% dither = [0,0.001:0.05:0.5].';
drift = [0.5*dVmean, dV_mean , 2*dVmean]/10;
initial_command = out.dm1.Vall(:,:,end);


% Pick the middle wavelength
si = mp.si_ref;%ceil(mp.Nsbp / 2);

contrast = zeros(n_iter,length(drift));


for i = 1:length(drift)
     % Initialize voltage vector
        
    V = zeros(mp.dm1.NactTotal,1);
    
    for j = 1:n_iter
       
        
        % set random command for active actuators?
        V(mp.dm1.act_ele) = V(mp.dm1.act_ele) + normrnd(0,drift(i),size(mp.dm1.act_ele));
        
        % Update mp.dm1 object
        mp.dm1.V = initial_command + reshape(V,[mp.dm1.Nact,mp.dm1.Nact]);
        
        % Take image with random command
        image = falco_get_sbp_image(mp,si); % this returns image in contrast units
        
        contrast(j,i) = mean(image(mp.Fend.score.mask));
        
        subplot(1,2,1)
        imagesc(log10(abs(image)))
        title(['image for $\sigma_{drift}$ = ',num2str(drift(i))],'Interpreter','latex')
        caxis([min(log10(abs(image(:)))), max(log10(abs(image(:))))])
        colorbar

        subplot(1,2,2)
        semilogy(1:j,contrast(1:j,i))
        xlim([1,50])
        ylim([0.1*min(contrast(:,i)),2*max(contrast(:,i))])
        xlabel('image number')
        ylabel('mean dark zone contrast')
        title([{'mean dark zone contrast'}, {'for $\sigma_{drift}$ = ',num2str(drift(i))}],'Interpreter','latex')

        drawnow
    end


end


%%


end_contrasts = contrast(:,end);
drift_legend = "drift = " + string(drift);


figure(232)
subplot(1,2,1)
semilogy(1:n_iter,contrast)
xlim([1,n_iter])
legend(drift_legend)
title([{'Mean dz contrast vs image number'},{'for varying dither'}])
xlabel('image number')
ylabel('mean dz contrast')


subplot(1,2,2)
plot(drift(2:end),end_contrasts)
xlabel('$\sigma_{dither} [V/\sqrt{iter}]$','Interpreter','latex')
ylabel('mean dz contrast')
title([{'mean dz contrast'},{'vs drift std dev'}])

%% Save
out_dither_test.contrast = contrast;
out_dither_test.mean_contrast = end_contrasts;
out_dither_test.drift = drift;

save(fullfile(dhmOutDir,['out_drift_test_',runTag,'.mat']),'out_dither_test')

