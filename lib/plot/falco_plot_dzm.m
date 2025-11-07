function hEplot = falco_plot_dzm(mp,ev)
        si = 1;
        fst = 20;
        Im = ev.imageArray(:,:,1,si);

        I_co_est = zeros(mp.Fend.Neta, mp.Fend.Nxi);
        I_co_est(mp.Fend.corr.mask) = abs(ev.Eest(:,si)).^2;

        I_inco_est = zeros(mp.Fend.Neta, mp.Fend.Nxi);
        I_inco_est(mp.Fend.corr.mask) = ev.IincoEst(:,si);

        %%-- Probed E-field plots
        hEplot = figure(222);
        set(hEplot,'units', 'inches', 'Position', [0 0 8 6])
        set(hEplot,'Color','w')

        h_psf = subplot(1,3,1); % Save the handle of the subplot
        % set(h_psf, 'OuterPosition', [0.01, 0.46, subplotbox])
        imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im),[-10 -3]); 
        axis xy equal tight; colorbar(h_psf); colormap(h_psf,parula);
        xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
        % ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
        set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
        title('$I$','Fontsize',fst,'Interpreter','LaTeX');

        h_est = subplot(1,3,2); % Save the handle of the subplot
        % set(h_est, 'OuterPosition', [0.01, 0.46, subplotbox])
        imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(I_co_est)); 
        axis xy equal tight; colorbar(h_est); colormap(h_est,parula);
        xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
        % ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
        set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
        title('$|\hat{E}_{co}|^2$','Fontsize',fst,'Interpreter','LaTeX');


        h_est_inco = subplot(1,3,3); % Save the handle of the subplot
        % set(h_est_inco, 'OuterPosition', [0.01, 0.46, subplotbox])
        imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(I_co_est)); 
        axis xy equal tight; colorbar(h_est_inco); colormap(h_est_inco,parula);
        xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
        % ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
        set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
        title('$\hat{I}_{inco}$','Fontsize',fst,'Interpreter','LaTeX');
        drawnow;

end