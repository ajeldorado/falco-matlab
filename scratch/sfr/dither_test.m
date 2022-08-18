%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test dithers to find contrast-dither relationship on IACT.
% Generate dark zone using preferred method and then run this script.
%
% sfr 12072022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dV_mean = mean(mean(abs(out.dm1.Vall(:,:,end) - out.dm1.Vall(:,:,end-1))));
p1 = 5.587e-6;
p2 = 0.14088e-6;
dx_mean = 2*p1*dV_mean + p2;
disp("mean displacement of DM actuators between last two iterations",num2str(dx_mean))
%%

figure(1)

% dither = [0,0.001:0.05:0.5].';
dither = [0, 0.5*dV_mean:dV_mean:dV_mean*10];
initial_command = out.dm1.Vall(:,:,end);


% Pick the middle wavelength
si = ceil(mp.Nsbp / 2);

contrast = zeros(50,size(dither,1));


for i = 1:1:size(dither,1)
    
    for j = 1:1:50
        % Initialize voltage vector
        V = zeros(mp.dm1.NactTotal,1);
        
        % set random command for active actuators?
        V(mp.dm1.act_ele) = normrnd(0,dither(i),size(mp.dm1.act_ele));
        
        % Update mp.dm1 object
        mp.dm1.V = initial_command + reshape(V,[mp.dm1.Nact,mp.dm1.Nact]);
        
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


%%


mean_contrasts = mean(contrast(:,2:end),1);
dither_legend = "dither = " + string(dither);


figure(2)
subplot(1,2,1)
semilogy([1:1:50],contrast)
xlim([1,50])
legend(dither_legend)
title([{'Mean dz contrast vs image number'},{'for varying dither'}])
xlabel('image number')
ylabel('mean dz contrast')


subplot(1,2,2)
plot(dither,mean_contrasts)
xlabel('$\sigma_{dither} [V/\sqrt{iter}]$','Interpreter','latex')
ylabel('mean dz contrast')
title([{'mean dz contrast'},{'vs dither std dev'}])

%% Save
out_dither_test.contrast = contrast;
out_dither_test.mean_contrast = mean_contrasts;
out_dither_test.dither = dither;

save('out_dither_test.mat',"out_dither_test")



