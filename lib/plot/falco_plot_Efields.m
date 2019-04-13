% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function handles = falco_plot_Efields(handles,mp,Itr,contrast_bandavg,Im,E)

IF3 = abs(E.F3).^2;
IF3 = IF3/max(IF3(:));

IP1 = abs(E.P1).^2;
IP3 = abs(E.P3).^2;
IP4 = abs(E.P4).^2;

if(mp.flagPlot)
    
    fig_size = [100 60 1500 900];

    if(Itr>1)
        delete(handles.tb1)
        delete(handles.tb2)
        delete(handles.tb3)
        delete(handles.tb4)
        figure(handles.master);
    else
        handles.master = figure('Color',[1 1 1],'Position',fig_size);
        set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    end

    h_text = subplot(2,3,1); % Save the handle of the subplot
    ax1=get(h_text,'position'); % Save the position as ax
    set(h_text,'position',ax1); % Manually setting this holds the position with colorbar 
    subplot(2,3,1); 
    axis off
    handles.tb1 = text(0.1,0.9,sprintf('%s: Iteration %d',mp.coro,Itr-1),'Fontsize',24);
    handles.tb2 = text(0.1,0.65,sprintf('%.1f%% BW @ %dnm',(100*mp.fracBW),round(mp.lambda0*1e9)),'Fontsize',24);
    handles.tb3 = text(0.1,0.35,sprintf('I_{norm} = %.2e',contrast_bandavg(Itr)),'Fontsize',24);
    handles.tb4 = text(0.1,0.05,sprintf('T_{core} =   %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',24);

    %--Middle top: F3
    h_psf = subplot(2,3,2); % Save the handle of the subplot
    ax2=get(h_psf,'position'); % Save the position as ax
    set(h_psf,'position',ax2); % Manually setting this holds the position with colorbar 
    subplot(2,3,2); 
    imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(IF3),[-5 0]); 
    axis xy equal tight; ch_psf=colorbar; colormap parula;
    xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
    ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    ylabel(ch_psf,'$log_{10}$(NI)','Fontsize',24,'Interpreter','LaTex');
    title('PSF at FPM','Fontsize',20,'Fontweight','Bold');

    %--Middle right: F4
    h_psf = subplot(2,3,3); % Save the handle of the subplot
    ax2=get(h_psf,'position'); % Save the position as ax
    set(h_psf,'position',ax2); % Manually setting this holds the position with colorbar 
    subplot(2,3,3); 
    imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im),[-10 -3]); 
    axis xy equal tight; ch_psf=colorbar; colormap parula;
    xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
    ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    ylabel(ch_psf,'$log_{10}$(NI)','Fontsize',24,'Interpreter','LaTex');
    title('Stellar PSF after Coronagraph','Fontsize',20,'Fontweight','Bold');

    %--Bottom Left: P1
    h_psf = subplot(2,3,4); % Save the handle of the subplot
    ax2=get(h_psf,'position'); % Save the position as ax
    set(h_psf,'position',ax2); % Manually setting this holds the position with colorbar 
    subplot(2,3,4); 
    imagesc((IP1)); 
    axis xy equal tight; colormap parula;   axis off; colorbar;
    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    title('Entrance Pupil','Fontsize',20,'Fontweight','Bold');
    
    %--Bottom Middle: P3
    h_psf = subplot(2,3,5); % Save the handle of the subplot
    ax2=get(h_psf,'position'); % Save the position as ax
    set(h_psf,'position',ax2); % Manually setting this holds the position with colorbar 
    subplot(2,3,5); 
    imagesc((IP3)); 
    axis xy equal tight; colormap parula;   axis off; colorbar;
    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    title('Pupil after DMs','Fontsize',20,'Fontweight','Bold');
    
    %--Bottom Right: P4
    h_psf = subplot(2,3,6); % Save the handle of the subplot
    ax2=get(h_psf,'position'); % Save the position as ax
    set(h_psf,'position',ax2); % Manually setting this holds the position with colorbar 
    subplot(2,3,6); 
    imagesc(log10(IP4),[-4 0]); 
    axis xy equal tight; colormap parula;   axis off; colorbar;
    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    title('Pupil at Lyot Plane','Fontsize',20,'Fontweight','Bold');
    
    drawnow;
   
end

end %--END OF FUNCTION