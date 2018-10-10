% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function handles = falco_plot_progress(handles,mp,Itr,contrast_bandavg,Im,DM1surf,DM2surf)


if(mp.flagPlot)
    fig_size = [100          60        1500         900];

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

    
    %--Need to specify slightly different coordinates for the SPHLC's FPM
    switch mp.coro
        case{'SPHLC'}
            mp.F3.compact.xisDL  = mp.F3.compact.in.xisDL;
            mp.F3.compact.etasDL = mp.F3.compact.in.etasDL;
    end
    
    
    switch mp.coro
        case{'Vortex','vortex','AVC','VC'}
            %--No subplot
        case{'EHLC'}
            h_amp = subplot(2,3,3);
             ax3=get(h_amp,'position'); % Save the position as ax
            set(h_amp,'position',ax3); % Manually setting this holds the position with colorbar 
            subplot(2,3,3); 
            imagesc(mp.dm8.compact.xisDL, mp.dm8.compact.xisDL, mp.DM8surf*1e9); axis xy equal tight; colormap gray; ch = colorbar;
            title('FPM Nickel Profile','Fontsize',16,'Fontweight','Bold');
            xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            ylabel(ch,'nm','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
        case{'HLC','APHLC','SPHLC'}
            h_amp = subplot(2,3,3);
             ax3=get(h_amp,'position'); % Save the position as ax
            set(h_amp,'position',ax3); % Manually setting this holds the position with colorbar 
            subplot(2,3,3); 
            imagesc(mp.F3.compact.xisDL,mp.F3.compact.xisDL,mp.DM8surf*1e9); axis xy equal tight; colormap gray; ch = colorbar;
            title('FPM Nickel Profile','Fontsize',16,'Fontweight','Bold');
            xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            ylabel(ch,'nm','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
        otherwise
            h_amp = subplot(2,3,3);
             ax3=get(h_amp,'position'); % Save the position as ax
            set(h_amp,'position',ax3); % Manually setting this holds the position with colorbar 
            subplot(2,3,3); 
            imagesc(mp.F3.compact.xisDL,mp.F3.compact.xisDL,mp.F3.compact.mask.amp); axis xy equal tight; ch = colorbar;
            title(sprintf('FPM Amplitude\n (%.2f%% ampl. transmission)',100*mp.FPMampFac),'Fontsize',16,'Fontweight','Bold');
            xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    end
    
    h_text = subplot(2,3,1); % Save the handle of the subplot
    ax1=get(h_text,'position'); % Save the position as ax
    set(h_text,'position',ax1); % Manually setting this holds the position with colorbar 
    subplot(2,3,1); 
    axis off
%     ht1 = text(0.05,0.9,sprintf('%s: Iteration %d',mp.coro,Itr-1),'Fontsize',24);
%     ht2 = text(0.05,0.5,sprintf('NI = %.2e',contrast_bandavg(Itr)),'Fontsize',24);
%     ht3 = text(0.05,0.1,sprintf('T_{core} = %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',24);
    handles.tb1 = text(0.1,0.9,sprintf('%s: Iteration %d',mp.coro,Itr-1),'Fontsize',24);
%     handles.tb1 = text(0.05,0.9,sprintf('%s, %d%% BW',mp.coro,round(100*mp.fracBW)),'Fontsize',24);
    handles.tb2 = text(0.1,0.65,sprintf('%.1f%% BW @ %dnm',(100*mp.fracBW),round(mp.lambda0*1e9)),'Fontsize',24);
    handles.tb3 = text(0.1,0.35,sprintf('I_{norm} = %.2e',contrast_bandavg(Itr)),'Fontsize',24);
    handles.tb4 = text(0.1,0.05,sprintf('T_{core} =   %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',24);

    h_psf = subplot(2,3,2); % Save the handle of the subplot
    ax2=get(h_psf,'position'); % Save the position as ax
    set(h_psf,'position',ax2); % Manually setting this holds the position with colorbar 
    subplot(2,3,2); 
    imagesc(mp.F4.xisDL,mp.F4.etasDL,log10(Im),[-10 -3]); 
    axis xy equal tight; ch_psf=colorbar; colormap parula;
    xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
    ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    ylabel(ch_psf,'$log_{10}$(NI)','Fontsize',24,'Interpreter','LaTex');
    title('Stellar PSF','Fontsize',20,'Fontweight','Bold');
%     title(sprintf('PSF, NI=%.2e',contrast_bandavg(Itr)),'Fontsize',20,'Fontweight','Bold');
%       title(sprintf('PSF at Iter=%03d,   NI=%.2e',Itr,contrast_bandavg(Itr)),'Fontsize',20,'Fontweight','Bold');


    h_dm1 = subplot(2,3,4); % Save the handle of the subplot
    ax4=get(h_dm1,'position'); % Save the position as ax
    set(h_dm1,'position',ax4); % Manually setting this holds the position with colorbar 
    subplot(2,3,4); 
    imagesc(1e9*DM1surf);  axis xy equal tight; axis off; ch1 = colorbar;
   ylabel(ch1,'nm','Fontsize',24,'Interpreter','LaTex');
%    title(sprintf('DM1 Surface at Iter=%03d',Itr),'Fontsize',20,'Fontweight','Bold');
   title(sprintf('DM1 Surface'),'Fontsize',20,'Fontweight','Bold');
   set(gca,'FontSize',20 ,'FontName','Times','FontWeight','Normal')

    h_dm2 = subplot(2,3,5); % Save the handle of the subplot
    ax5=get(h_dm2,'position'); % Save the position as ax
    set(h_dm2,'position',ax5); % Manually setting this holds the position with colorbar 
    subplot(2,3,5); 
    imagesc(1e9*DM2surf);  axis xy equal tight; axis off; ch1 = colorbar;
    ylabel(ch1,'nm','Fontsize',24,'Interpreter','LaTex');
   title(sprintf('DM2 Surface'),'Fontsize',20,'Fontweight','Bold');
    set(gca,'FontSize',20 ,'FontName','Times','FontWeight','Normal')

    switch mp.coro
%         case{'Vortex'}
%             %--No subplot
        case{'EHLC'}
            
           
            h_dm9 = subplot(2,3,6);
            ax6=get(h_dm9,'position'); % Save the position as ax
            set(h_dm9,'position',ax6); % Manually setting this holds the position with colorbar 
            subplot(2,3,6); 
            imagesc(mp.dm9.compact.xisDL, mp.dm9.compact.xisDL, mp.DM9surf*1e9); axis xy equal tight; ch = colorbar;
            title('FPM PMGI Profile','Fontsize',16,'Fontweight','Bold');
            xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            ylabel(ch,'nm','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
        
        case{'HLC','APHLC','SPHLC'}
            
            if( max(mp.DM9surf(:))==0 )
                cutoffFloor = max(mp.DM9surf(:)) -1e-9;
            else
            
                DM9surfNZ = mp.DM9surf(mp.DM9surf~=0);
                DM9surfNZsort = sort(DM9surfNZ(:));
                cutoffFrac = 0.02;
                cutoffInd = ceil(cutoffFrac*length(DM9surfNZsort));
                if(cutoffInd<1)
                    cutoffInd = 1;
                elseif(cutoffInd>length(DM9surfNZsort))
                    cutoffInd = length(DM9surfNZsort);
                end
                cutoffFloor = DM9surfNZsort( cutoffInd );
                if( cutoffFloor >= max(mp.DM9surf(:)) )
                    cutoffFloor = 0.9*max(mp.DM9surf(:));
                end
                
            end
           
            h_dm9 = subplot(2,3,6);
            ax6=get(h_dm9,'position'); % Save the position as ax
            set(h_dm9,'position',ax6); % Manually setting this holds the position with colorbar 
            subplot(2,3,6); 
            imagesc(mp.F3.compact.xisDL,mp.F3.compact.xisDL,mp.DM9surf*1e9,1e9*[cutoffFloor,max(mp.DM9surf(:))]); axis xy equal tight; ch = colorbar;
            title('FPM PMGI Profile','Fontsize',16,'Fontweight','Bold');
            xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            ylabel(ch,'nm','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
        
        otherwise %--No subplot

    end
    
    
    
    drawnow;
   
end



end %--END OF FUNCTION
