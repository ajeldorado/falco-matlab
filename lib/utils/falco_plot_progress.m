% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function handles = falco_plot_progress(handles,mp,Itr,contrast_bandavg,ImBandAvg_array,DM1S_array,DM2S_array)


if(mp.flagPlot)
    fig_size = [100          100        1200         800];

    if(Itr>1)
        delete(handles.tb1)
        delete(handles.tb2)
        delete(handles.tb3)
        delete(handles.tb4)
        figure(handles.master);
    else
%         handles.master = figure(123);
%         set(gcf,'Color',[1 1 1],'Position',fig_size);
        handles.master = figure('Color',[1 1 1],'Position',fig_size);
        set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    end

    
    h_text = subplot(2,2,1); % Save the handle of the subplot
    ax1=get(h_text,'position'); % Save the position as ax
    set(h_text,'position',ax1); % Manually setting this holds the position with colorbar 
    subplot(2,2,1); 
    axis off
%     ht1 = text(0.10,0.9,sprintf('%s: Iteration %d',mp.coro,Itr-1),'Fontsize',24);
%     ht2 = text(0.10,0.5,sprintf('NI = %.2e',contrast_bandavg(Itr)),'Fontsize',24);
%     ht3 = text(0.10,0.1,sprintf('T_{core} = %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',24);
    handles.tb1 = text(0.10,0.9,sprintf('%s: Iteration %d',mp.coro,Itr-1),'Fontsize',24);
%     handles.tb1 = text(0.10,0.9,sprintf('%s, %d%% BW',mp.coro,round(100*mp.fracBW)),'Fontsize',24);
    handles.tb2 = text(0.10,0.65,sprintf('%d%% BW @ %dnm',round(100*mp.fracBW),round(mp.lambda0*1e9)),'Fontsize',24);
    handles.tb3 = text(0.10,0.35,sprintf('I_{norm} = %.2e',contrast_bandavg(Itr)),'Fontsize',24);
    handles.tb4 = text(0.10,0.05,sprintf('T_{core} =   %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',24);

    h_psf = subplot(2,2,2); % Save the handle of the subplot
    ax2=get(h_psf,'position'); % Save the position as ax
    set(h_psf,'position',ax2); % Manually setting this holds the position with colorbar 
    subplot(2,2,2); 
    imagesc(mp.F4.full.xisDL,mp.F4.full.etasDL,log10(ImBandAvg_array(:,:,Itr)),[-10 -3]); 
    axis xy equal tight; ch_psf=colorbar; colormap parula;
    xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
    ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    ylabel(ch_psf,'$log_{10}$(NI)','Fontsize',20,'Interpreter','LaTex');
    title('Stellar PSF','Fontsize',20,'Fontweight','Bold');
%     title(sprintf('PSF, NI=%.2e',contrast_bandavg(Itr)),'Fontsize',20,'Fontweight','Bold');
%       title(sprintf('PSF at Iter=%03d,   NI=%.2e',Itr,contrast_bandavg(Itr)),'Fontsize',20,'Fontweight','Bold');


    h_dm1 = subplot(2,2,3); % Save the handle of the subplot
    ax4=get(h_dm1,'position'); % Save the position as ax
    set(h_dm1,'position',ax4); % Manually setting this holds the position with colorbar 
    subplot(2,2,3); 
    imagesc(DM1S_array(:,:,Itr)*1e9);  axis xy equal tight; axis off; ch1 = colorbar;
   ylabel(ch1,'nm','Fontsize',20,'Interpreter','LaTex');
%    title(sprintf('DM1 Surface at Iter=%03d',Itr),'Fontsize',20,'Fontweight','Bold');
   title(sprintf('DM1 Surface'),'Fontsize',20,'Fontweight','Bold');
   set(gca,'FontSize',20 ,'FontName','Times','FontWeight','Normal')

    h_dm2 = subplot(2,2,4); % Save the handle of the subplot
    ax5=get(h_dm2,'position'); % Save the position as ax
    set(h_dm2,'position',ax5); % Manually setting this holds the position with colorbar 
    subplot(2,2,4); 
    imagesc(DM2S_array(:,:,Itr)*1e9);  axis xy equal tight; axis off; ch1 = colorbar;
    ylabel(ch1,'nm','Fontsize',20,'Interpreter','LaTex');
   title(sprintf('DM2 Surface'),'Fontsize',20,'Fontweight','Bold');
    set(gca,'FontSize',20 ,'FontName','Times','FontWeight','Normal')


   
end



end %--END OF FUNCTION