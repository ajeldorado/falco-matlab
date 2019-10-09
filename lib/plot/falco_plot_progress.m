% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function handles = falco_plot_progress(handles,mp,Itr,InormHist,Im,DM1surf,DM2surf,ImSimOffaxis)

Im(Im<0) = 0; %--Prevent the log10(Im) plot from getting complex values.

if(mp.flagPlot)
    
    %--Figure is nominally a square. Need extra column when mp.coro='HLC'
    switch upper(mp.coro)
        case{'HLC'}
            fw = 1000;
            fh = 700;
            fig_size = [10 10 fw fh];
        otherwise
            fig_size = [10 10 800 800];
    end

    if(Itr>1)
        delete(handles.textbox1); %--Don't want to overlay the top label
        try
            figure(handles.master);
        catch
            handles.master = figure('Color',[1 1 1],'Position',fig_size);
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
        end
    else
        handles.master = figure('Color',[1 1 1],'Position',fig_size);
        set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
    end

    switch upper(mp.coro)
        case{'HLC'}
            
            subplotbox = [1, fw/fh]*0.32;
            fst = 20; %--Font size for titles in the subplots

            %--Top label
            handles.textbox1 = annotation('textbox', [0.25 0.89 0.5 0.1], ...
                'String', sprintf('%s: Iteration %d',mp.coro,Itr-1),'Fontsize',32,...
                'HorizontalAlignment','center','LineStyle','none','Interpreter','latex');

            %--Stellar PSF
            h_psf = subplot(2,3,1); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.01, 0.46, subplotbox])
            imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im),[-10 -3]); 
            axis xy equal tight; colorbar(h_psf); colormap(h_psf,parula);
            % xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
            %ylabel(ch_psf,'$log_{10}$(NI)','Fontsize',24,'Interpreter','LaTex');
            title(sprintf('Stellar PSF: NI = %.2e',InormHist(Itr)),'Fontsize',fst);%,'Fontweight','Bold');

            %--Off-axis PSF
            h_offaxis = subplot(2,3,4); % Save the handle of the subplot
            set(h_offaxis, 'OuterPosition', [0.01, 0.02, subplotbox])
            %--Crop the plot to be around the off-axis PSF
            [XIS,ETAS] = meshgrid(mp.Fend.eval.xisDL,mp.Fend.eval.etasDL);
            [~,ind_max] = max(ImSimOffaxis(:));
            xi_max = XIS(ind_max);
            eta_max = ETAS(ind_max);
            buffer = 4; %--lambda/D
            xi_range = buffer*[-1 1] + xi_max;
            eta_range = buffer*[-1 1] + eta_max;
            imagesc(mp.Fend.eval.xisDL,mp.Fend.eval.etasDL,ImSimOffaxis/(max(ImSimOffaxis(:))),[0 1]); 
            axis xy equal tight; colorbar; colormap(h_offaxis,'parula'); %colormap(h_offaxis,brewermap([],'Blues'));
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
            xlim(xi_range);
            ylim(eta_range);
            switch lower(mp.thput_metric)
                case{'hmi'} %--Absolute energy within half-max isophote(s)
                    title(sprintf('Off-axis PSF: T_{HM} =   %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',fst);
                case{'ee','e.e.'} %--Absolute energy encircled within a given radius
                    title(sprintf('Off-axis PSF: T_{EE} =   %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',fst);
            end

            h_dm1 = subplot(2,3,2); % Save the handle of the subplot
            set(h_dm1, 'OuterPosition', [0.35, 0.46, subplotbox])
            imagesc(1e9*DM1surf);  axis xy equal tight; axis off; colorbar; colormap(h_dm1,'parula');
            %ylabel(ch1,'nm','Fontsize',24,'Interpreter','LaTex');
            title(sprintf('DM1 Surface (nm)'),'Fontsize',fst) %,'Fontweight','Bold');
            set(gca,'FontSize',20 ,'FontName','Times','FontWeight','Normal')

            h_dm2 = subplot(2,3,5); % Save the handle of the subplot
            set(h_dm2, 'OuterPosition', [0.35, 0.02, subplotbox])
            imagesc(1e9*DM2surf);  axis xy equal tight; axis off; colorbar; colormap(h_dm2,'parula');
            %ylabel(ch1,'nm','Fontsize',24,'Interpreter','LaTex');
            title(sprintf('DM2 Surface (nm)'),'Fontsize',fst) %,'Fontweight','Bold');
            set(gca,'FontSize',20 ,'FontName','Times','FontWeight','Normal')

            switch lower(mp.layout)
                case 'fourier'
                    h_dm8 = subplot(2,3,3);
                    set(h_dm8, 'OuterPosition', [0.67, 0.46, subplotbox])
                    imagesc(mp.F3.compact.xisDL,mp.F3.compact.xisDL,mp.DM8surf*1e9); axis xy equal tight; colormap parula; colorbar;
                    title('FPM Nickel (nm)','Fontsize',16,'Fontweight','Bold');
%                     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
                    ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
%                     ylabel(ch,'nm','FontSize',16,'Interpreter','LaTeX');
                    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
                    
                    %--Don't plot down to zero in order to see the features
                    % on the dielectric profile.
                    if( max(mp.DM9surf(:))==0 )
                        cutoffFloor = max(mp.DM9surf(:)) -1e-9;
                    else

                        DM9surfNZ = mp.DM9surf(mp.DM9surf~=0);
                        DM9surfNZsort = sort(DM9surfNZ(:));
                        cutoffFrac = 0.04;
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
                    set(h_dm9, 'OuterPosition', [0.67, 0.02, subplotbox])
                    imagesc(mp.F3.compact.xisDL,mp.F3.compact.xisDL,mp.DM9surf*1e9,1e9*[cutoffFloor,max(mp.DM9surf(:))]); axis xy equal tight; colormap parula; colorbar;
                    title('FPM PMGI (nm)','Fontsize',16,'Fontweight','Bold');
        %             xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
                    ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
        %             ylabel(ch,'nm','FontSize',16,'Interpreter','LaTeX');
                    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
            end

        otherwise
            
            handles.textbox1 = annotation('textbox', [0.25 0.86 0.5 0.1], ...
                'String', sprintf('%s: Iteration %d',mp.coro,Itr-1),'Fontsize',32,...
                'HorizontalAlignment','center','LineStyle','none','Interpreter','latex');

            fst = 24; %--Font size for titles in the subplots

            h_psf = subplot(2,2,1); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.46, [1 1]*0.45])
            imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im),[-10 -3]); 
            axis xy equal tight; colorbar(h_psf); colormap(h_psf,parula);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
            %ylabel(ch_psf,'$log_{10}$(NI)','Fontsize',24,'Interpreter','LaTex');
            title(sprintf('Stellar PSF: NI = %.2e',InormHist(Itr)),'Fontsize',fst);%,'Fontweight','Bold');

            h_offaxis = subplot(2,2,3); % Save the handle of the subplot
            set(h_offaxis, 'OuterPosition', [0.05, 0.02, [1 1]*0.45])
            %--Crop the plot to be around the off-axis PSF
            [XIS,ETAS] = meshgrid(mp.Fend.eval.xisDL,mp.Fend.eval.etasDL);
            [~,ind_max] = max(ImSimOffaxis(:));
            xi_max = XIS(ind_max);
            eta_max = ETAS(ind_max);
            buffer = 4; %--lambda/D
            xi_range = buffer*[-1 1] + xi_max;
            eta_range = buffer*[-1 1] + eta_max;
            imagesc(mp.Fend.eval.xisDL,mp.Fend.eval.etasDL,ImSimOffaxis/(max(ImSimOffaxis(:))),[0 1]); 
            axis xy equal tight; colorbar; colormap(h_offaxis,'parula'); %colormap(h_offaxis,brewermap([],'Blues'));
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
            set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
            xlim(xi_range);
            ylim(eta_range);
            switch lower(mp.thput_metric)
                case{'hmi'} %--Absolute energy within half-max isophote(s)
                    title(sprintf('Off-axis PSF: T_{HM} =   %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',fst);
                case{'ee','e.e.'} %--Absolute energy encircled within a given radius
                    title(sprintf('Off-axis PSF: T_{EE} =   %.2f%%',100*mp.thput_vec(Itr)),'Fontsize',fst);
            end

            h_dm1 = subplot(2,2,2); % Save the handle of the subplot
            set(h_dm1, 'OuterPosition', [0.55, 0.46, [1 1]*0.45])
            imagesc(1e9*DM1surf);  axis xy equal tight; axis off; colorbar; colormap(h_dm1,'parula');
            %ylabel(ch1,'nm','Fontsize',24,'Interpreter','LaTex');
            title(sprintf('DM1 Surface (nm)'),'Fontsize',fst) %,'Fontweight','Bold');
            set(gca,'FontSize',20 ,'FontName','Times','FontWeight','Normal')

            h_dm2 = subplot(2,2,4); % Save the handle of the subplot
            set(h_dm2, 'OuterPosition', [0.55, 0.02, [1 1]*0.45])
            imagesc(1e9*DM2surf);  axis xy equal tight; axis off; colorbar; colormap(h_dm2,'parula');
            %ylabel(ch1,'nm','Fontsize',24,'Interpreter','LaTex');
            title(sprintf('DM2 Surface (nm)'),'Fontsize',fst) %,'Fontweight','Bold');
            set(gca,'FontSize',20 ,'FontName','Times','FontWeight','Normal')            
            
    end   

    drawnow;
   
end

end %--END OF FUNCTION