% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% handles = falco_plot_progress_omc(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf)

function handles = falco_plot_progress_omc(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf)

if mp.flagSim
    % handles = falco_plot_progress_omc_model(handles, mp, Itr, InormHist_tb, Im_tb, DM1surf, DM2surf);
    % handles = falco_plot_progress_omc_model(handles, mp, Itr, Inorm, Im_tb, DM1surf, DM2surf);
    % return
    tb = [];
else
    tb = mp.tb;
end

% only difference between testbed and model is tb writes some tb stuff to
% fits header
if isempty(tb)
    fFitsWrite = @(tb, im, fn) fitswrite(im, fn);
else
    fFitsWrite = @(tb, im, fn) sciCam_fitswrite(tb, im, fn);
end

subplot = @(m,n,p) subtightplot(m,n,p,[0.025 0.025],[0.1 0.1],[0.1 0.1]);

Icbmin = -10;
Icbmax = -4;

Im = Im_tb.Im;
Im4plot = Im;
Im4plot(Im4plot<0) = 0; %--Prevent the log10(Im) plot from getting complex values.

if(mp.flagPlot)

    if(Itr>1)
%         delete(handles.tb1)
%         delete(handles.tb2)
%         delete(handles.tb3)
%         delete(handles.tb4)
%         delete(handles.tb5)
        try
            figure(handles.master);
        catch
            handles.master = figure('Color','w');
            set(handles.master,'units', 'inches', 'Position', [0 0 12 8])
        end
    else
        handles.master = figure('Color','w');
        set(handles.master,'units', 'inches', 'Position', [0 0 12 8])
    end

%     subplot(2,3,1); % Save the handle of the subplot
%     axis off
%     handles.tb1 = text(0.1,0.8,sprintf('%s: Iteration %d',mp.coro,Itr-1));
%     handles.tb2 = text(0.1,0.7,sprintf('%.1f%% BW @ %dnm',(100*mp.fracBW),round(mp.lambda0*1e9)));
%     handles.tb3 = text(0.1,0.6,sprintf('I_{norm} = %.2e',Inorm.total(Itr)));
%     handles.tb4 = text(0.1,0.5,sprintf('I_{mod,prev} = %.2e',Imod));
%     switch lower(mp.thput_metric)
%         case{'hmi'} %--Absolute energy within half-max isophote(s)
%             handles.tb5 = text(0.1,0.4,sprintf('T_{half-max} =   %.2f%%',100*mp.thput_vec(Itr)));
%         case{'ee','e.e.'} %--Absolute energy encircled within a given radius
%             handles.tb5 = text(0.1,0.4,sprintf('T_{E.E.} =   %.2f%%',100*mp.thput_vec(Itr)));
%     end

    subplot(2,3,1); 
    imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im4plot),[Icbmin Icbmax]); 
    axis xy equal tight; 
    colorbar; 
    colormap(gca,parula);
    try
        title(['it = ',num2str(Itr-1),', Inorm = ',num2str(Inorm.total(Itr),2)]);
    catch % sometimes the Inorm total doesnt have a new value, like the last plot update
        title(['it = ',num2str(Itr-1),', Inorm = ',num2str(Inorm.total(Itr-1),2)]);
    end
    
    %% add rectangular correction / scoring regions
    if(exist(mp.Fend.shape)); 
        if strcmp(mp.Fend.shape,'square')
            hold on; rectangle('Position', [mp.Fend.xiOffset - mp.Fend.corr.Rout, mp.Fend.etaOffset - mp.Fend.corr.Rout, mp.Fend.corr.Rout*2, mp.Fend.corr.Rout*2],'EdgeColor','r');
            hold on; rectangle('Position', [mp.Fend.xiOffset - mp.Fend.score.Rout, mp.Fend.etaOffset - mp.Fend.score.Rout, mp.Fend.score.Rout*2, mp.Fend.score.Rout*2],'EdgeColor','y');
        end
    end
        
    try
        axis(mp.Fend.dzAxis)
    catch
        axis xy equal tight; 
    end
    

	subplot(2,3,2); 
    imagesc(1e9*DM1surf);  axis xy equal tight; axis off;
    colorbar;
    colormap(gca,gray);
    title('DM1 Surface (nm)');

	subplot(2,3,3); 
    imagesc(1e9*DM2surf);  axis xy equal tight; axis off;
    colorbar;
    colormap(gca,gray);
    title('DM2 Surface (nm)');

    subplot(2,3,4);
    semilogy(0:length(Inorm.total)-1,Inorm.total,'-o');hold on;
    semilogy(0:Itr-1,mean(Inorm.mod,2),'-o');
    semilogy(0:Itr-1,mean(Inorm.unmod,2),'--o');
    hold off;
    xlim([0 length(Inorm.total)])
    xlabel('Iteration')
%     ylabel('Norm. I');
    legend('Total','Modulated','Unmodulated');
	title('Mean Normalized Intensity')
    grid on;axis square;
% 	hcbdummy = colorbar;set(hcbdummy,'visible','off');
    
	subplot(2,3,5)
    cmap = jet(mp.Nsbp+1);
    cmap = cmap ./ sum(cmap,2);% make the jet cmap darker

    legstr = {};
    for istar = 1:mp.star.count
        for isb = 1:mp.Nsbp
            si = isb + (istar-1)*mp.Nsbp;
            if(si==mp.si_ref)
                linecolor=[0 0 0];
            elseif(si==mp.Nsbp)
                linecolor=cmap(end,:);% force last band to red
            else
                linecolor=cmap(si,:);
            end
            legstr{si} = ['sbnd ' num2str(isb) '; star ' num2str(istar)];
            semilogy(0:Itr-1,Inorm.mod(:,si),'-o','Color',linecolor); hold on;
            %hl2(si)=semilogy(0:Itr-2,Inorm.unmod(:,si),'--o','Color',linecolor);
        end % subband
    end % star
    legend(legstr{:});
    hold off;
    xlim([0 length(Inorm.total)])
    xlabel('Iteration')
%     ylabel('Norm. I');
    %if(Itr>2); legend([hl1(mp.si_ref), hl2(mp.si_ref)],'Modulated','Unmodulated');end
	title('Mean Mod Intensity')
    grid on;axis square;
% 	hcbdummy = colorbar;set(hcbdummy,'visible','off');

	subplot(2,3,6)

    if(mp.Nsbp>1)
        semilogy(mp.sbp_centers*1e9,Inorm.mod(end,:),'-o');
    else
        semilogy(mp.sbp_centers*1e9,Inorm.mod(end),'-o');
    end

%     hold off;
    xlabel('Wavelength (nm)')
%     legend('Mean Total','Modulated','location','best');
	title('Mean Mod Intensity')
    grid on;axis square;
% 	hcbdummy = colorbar;set(hcbdummy,'visible','off');

% 	subplot(2,3,5); % Save the handle of the subplot
%     imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(abs(Im_tb.E(:,:,si_ref)).^2),[Icbmin Icbmax]); 
%     axis xy equal tight;
%     colorbar;
%     colormap(gca,parula)  
% %     xlabel('\lambda_0/D'); 
% %     ylabel('\lambda_0/D');
%     title('Modulated (previous)');
%     
% 	subplot(2,3,6); % Save the handle of the subplot
%     imagesc(mp.Fend.xisDL,mp.Fend.etasDL,angle(Im_tb.E(:,:,si_ref)),[-pi pi]); 
%     axis xy equal tight; 
%     colorbar; 
%     colormap(gca,hsv);
% %     xlabel('\lambda_0/D'); 
% %     ylabel('\lambda_0/D');
%     title('Phase (previous)');
%     
   FigureTitle(['Trial ' num2str(mp.TrialNum, '%04d')]);
   drawnow;


    if(Itr>1)
        %%-- Probed E-field plots
        hEplot = figure(98);
        set(hEplot,'units', 'inches', 'Position', [0 0 4*mp.Nsbp 6])
        set(hEplot,'Color','w')

        for si = 1:mp.Nsbp
            subplot(2,mp.Nsbp,si); % Save the handle of the subplot
            imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(abs(Im_tb.E(:,:,si)).^2),[Icbmin Icbmax]); 
            axis xy equal tight;
            colorbar;
            colormap(gca,parula)  
            try
                axis(mp.Fend.dzAxis)
            catch
                axis xy equal tight;
            end
            title(['|E_{tb}|, band ' num2str(si)])
    
            subplot(2,mp.Nsbp,si+mp.Nsbp); % Save the handle of the subplot
            imagesc(mp.Fend.xisDL,mp.Fend.etasDL,angle(Im_tb.E(:,:,si)),[-pi pi]); 
            axis xy equal tight; 
            colorbar; 
            colormap(gca,hsv);
            title(['Arg(E_{tb}), band ' num2str(si)])
    
            try
                 axis(mp.Fend.dzAxis)
            catch
                 axis xy equal tight; 
            end
        end
        drawnow;
    end

	%%-- TO DO: 
    % - DM stroke usage (rms,ptv) plots 
    % - Throughput plots 
    % - 

end % if plot


%%-- Save data
if isempty(tb)
    % mp.OUT_DATA_DIR = fullfile(getenv("DATA_ROOT"), ['falco_testbed_run' num2str(mp.SeriesNum)], 'data', mp.runLabel);
    out_dir = fullfile(mp.OUT_DATA_DIR, mp.runLabel);
else
    out_dir = fullfile(tb.info.OUT_DATA_DIR, mp.runLabel);
end
% Directory to save dat
if(~exist(out_dir, 'dir'))
    mkdir(out_dir);
end

fFitsWrite(tb,Im,fullfile(out_dir,['normI_it',num2str(Itr-1),'.fits']));

if(any(mp.dm_ind==1) && Itr==1)
    fFitsWrite(tb,mp.dm1.biasMap,fullfile(out_dir,'dm1_Vbias.fits'));
end
if(any(mp.dm_ind==2) && Itr==1)
    fFitsWrite(tb,mp.dm2.biasMap,fullfile(out_dir,'dm2_Vbias.fits'));
end

if(any(mp.dm_ind==1))
    fFitsWrite(tb,mp.dm1.V,fullfile(out_dir,['dm1_V_it',num2str(Itr-1),'.fits']));
    fFitsWrite(tb,DM1surf,fullfile(out_dir,['dm1_model_it',num2str(Itr-1),'.fits']));
end
if(any(mp.dm_ind==2))
    fFitsWrite(tb,mp.dm2.V,fullfile(out_dir,['dm2_V_it',num2str(Itr-1),'.fits']));
    fFitsWrite(tb,DM2surf,fullfile(out_dir,['dm2_model_it',num2str(Itr-1),'.fits']));
end


% Im_tb.E should be a cube with all modes but sciCam_fitswrite does not
% write it as a cube
%need to also export the unprobed images
% for iStar = 1:mp.compact.star.count
%     for si = 1:mp.Nsbp
%         modeIndex = (iStar-1)*mp.Nsbp + si;
%                     
%         %reference
%         %tmp = zeros(size(Im));
%         %tmp(mp.Fend.corr.maskBool) = ev.Eest(:, modeIndex);
%         %Im_tb.E(:, :, modeIndex) = tmp; % modulated component
%         %tmp = zeros(size(Im));
%         %tmp(mp.Fend.corr.maskBool) = ev.IincoEst(:, si);
%         %Im_tb.Iinco(:, :, modeIndex) = tmp; % unmodulated component
%         %out.InormHist_tb.mod(Itr, modeIndex) = mean(abs(ev.Eest(:, modeIndex)).^2);
%         %out.InormHist_tb.unmod(Itr, modeIndex) = mean(ev.IincoEst(:, modeIndex)); 
%         
%         
%         thisIunprobed = Im_tb.ev.I0{modeIndex};
%         
%         %tmp = zeros(size(Im));
%         %tmp = squeeze(Im_tb.E(:,:,modeIndex));
%         %thisEsens = zeros(size(Im));
%         thisEsens = squeeze(Im_tb.E(:,:,modeIndex));
%         
%         %tmp = zeros(size(Im));
%         %tmp = squeeze(Im_tb.Iinco(:,:,modeIndex));
%         %thisIinco = zeros(size(Im));
%         %thisIinco(mp.Fend.corr.maskBool) = tmp;
%         thisIinco = squeeze(Im_tb.Iinco(:,:,modeIndex));
%         
%             
%         %new
%         fFitsWrite(tb,abs(thisEsens).^2,fullfile(out_dir,['normI_Esens_it',num2str(Itr-1), '_mode', num2str(modeIndex),'.fits']));
%         fFitsWrite(tb,angle(thisEsens),fullfile(out_dir,['phz_Esens_it',num2str(Itr-1), '_mode', num2str(modeIndex),'.fits']));
%         fFitsWrite(tb,thisIinco,fullfile(out_dir,['normI_inco_it',num2str(Itr-1), '_mode', num2str(modeIndex),'.fits']));
%         fFitsWrite(tb,thisIunprobed,fullfile(out_dir,['normI_unprobed_it',num2str(Itr-1),'_mode', num2str(modeIndex),'.fits']));
%      end
% end

%add unprobed fits files
for iStar = 1:mp.compact.star.count
    for si = 1:mp.Nsbp
        modeIndex = (iStar-1)*mp.Nsbp + si;
        thisIunprobed = Im_tb.ev.I0{modeIndex};
    end
end
fFitsWrite(tb,thisIunprobed,fullfile(out_dir,['normI_unprobed_it',num2str(Itr-1),'_mode', num2str(modeIndex),'.fits']));


%old
fFitsWrite(tb,abs(Im_tb.E).^2,fullfile(out_dir,['normI_Esens_it',num2str(Itr-1),'.fits']));
fFitsWrite(tb,angle(Im_tb.E),fullfile(out_dir,['phz_Esens_it',num2str(Itr-1),'.fits']));
fFitsWrite(tb,Im_tb.Iinco,fullfile(out_dir,['normI_inco_it',num2str(Itr-1),'.fits']));
 

% all the data is being saved in the .mat structure
if(~strcmpi(mp.estimator,'perfect'))
    ev = Im_tb.ev;
    save(fullfile(out_dir,['probing_data_',num2str(Itr-1),'.mat']),'ev');
end

% % Update the diary 
% diary off; diary(mp.diaryfile)

end %--END OF FUNCTION

function han_out = FigureTitle(stitle, varargin)
% han = FigureTitle(stitle, varargin)
%
% add an annotation text at the top of the figure, 
% useful for adding a single main title to a figure with subplots
%
% varargin can be property, value pairs sent to the annotation handle
% copied from D. Marx matlab toolbox

han = annotation('textbox', [0.5 0.8 0.2 0.2], 'String', stitle, ...
    'FitBoxToText', 'on', 'LineStyle', 'none', ...
    'FontSize', 24, 'Color', 'r', 'FontWeight', 'bold');
set(han,'HorizontalAlignment','center')
% center horizontally
ppp = get(han,'Position');
set(han,'Position',[0.5 - 0.5*ppp(3) ppp(2:end)])
% so it can be found and deleted later
set(get(han,'parent'),'HandleVisibility','on')

if ~isempty(varargin),
    set(han, varargin{:})
end

if nargout > 0,
    han_out = han;
end

end % FigureTitle