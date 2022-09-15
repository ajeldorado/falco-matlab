% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function handles = falco_plot_progress_dst(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf)

tb = mp.tb;

if(Itr==10 || Itr==40 && ~mp.flagSim)
    % Clear the dark 
    disp('Clearing the dark ...');
    sbp_texp = tb.info.sbp_texp(mp.si_ref);
    [~,flnm] = sciCam_loadDark(tb,sbp_texp);
    delete(flnm);
else
    disp('Keeping dark ...');
end


subplot = @(m,n,p) subtightplot(m,n,p,[0.025 0.025],[0.1 0.1],[0.1 0.1]);

Icbmin = -10;
Icbmax = -4;

Im = Im_tb.Im;
Im4plot = Im;
Im4plot(Im4plot<0) = 0; %--Prevent the log10(Im) plot from getting complex values.

% if(Itr>1 && ~strcmpi(mp.estimator,'perfect') )
%     Imod = mean(Inorm.mod(Itr-1,:));
% else
%     Imod = NaN;
% end



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
    if strcmpi(mp.estimator,'ekf_maintenance') %&& any(mp.est.itr_ol==Itr) == true
        % TODO: this fails when it isnt an OL iteration for some reason
        semilogy(0:Itr-1,mean(Im_tb.ev.IOLScoreHist(1:Itr,:),2),'-p')
    end
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

    for si = 1:mp.Nsbp
        if(si==mp.si_ref)
            linecolor=[0 0 0];
        elseif(si==mp.Nsbp)
            linecolor=cmap(end,:);% force last band to red
        else
            linecolor=cmap(si,:);
        end
        semilogy(0:Itr-1,Inorm.mod(:,si),'-o','Color',linecolor); hold on;
        %hl2(si)=semilogy(0:Itr-2,Inorm.unmod(:,si),'--o','Color',linecolor);
    end

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
        %     xlabel('\lambda_0/D'); 
        %     ylabel('\lambda_0/D');
    %         title(['Modulated (band ',num2str(si),')']);

            subplot(2,mp.Nsbp,si+mp.Nsbp); % Save the handle of the subplot
            imagesc(mp.Fend.xisDL,mp.Fend.etasDL,angle(Im_tb.E(:,:,si)),[-pi pi]); 
            axis xy equal tight; 
            colorbar; 
            colormap(gca,hsv);
        %     xlabel('\lambda_0/D'); 
        %     ylabel('\lambda_0/D');
    %         title(['Phase (band ',num2str(si),')']);
        end
        drawnow;
    end

	%%-- TO DO: 
    % - DM stroke usage (rms,ptv) plots 
    % - Throughput plots 
    % - 

end


%%-- Save data

out_dir = fullfile(tb.info.OUT_DATA_DIR,mp.runLabel);
% Directory to save dat
if(~exist(out_dir, 'dir'))
    mkdir(out_dir);
end

% Make it clear that the file is simulated if mp.flagSim = true
if(mp.flagSim)
    tag = '_SIM';
else
    tag = '';
end

fitswrite_tb(mp,tb,Im,fullfile(out_dir,['normI_it',num2str(Itr-1),tag,'.fits']));

if(any(mp.dm_ind==1) && Itr==1)
    fitswrite_tb(mp,tb,mp.dm1.biasMap,fullfile(out_dir,'dm1_Vbias.fits'));
end
if(any(mp.dm_ind==2) && Itr==1)
    fitswrite_tb(mp,tb,mp.dm2.biasMap,fullfile(out_dir,'dm2_Vbias.fits'));
end

if(any(mp.dm_ind==1))
    fitswrite_tb(mp,tb,mp.dm1.V,fullfile(out_dir,['dm1_V_it',num2str(Itr-1),tag,'.fits']));
    fitswrite_tb(mp,tb,DM1surf,fullfile(out_dir,['dm1_model_it',num2str(Itr-1),tag,'.fits']));
end
if(any(mp.dm_ind==2))
    fitswrite_tb(mp,tb,mp.dm2.V,fullfile(out_dir,['dm2_V_it',num2str(Itr-1),tag,'.fits']));
    fitswrite_tb(mp,tb,DM2surf,fullfile(out_dir,['dm2_model_it',num2str(Itr-1),tag,'.fits']));
end


fitswrite_tb(mp,tb,abs(Im_tb.E).^2,fullfile(out_dir,['normI_Esens_it',num2str(Itr-1),tag,'.fits']));
fitswrite_tb(mp,tb,angle(Im_tb.E),fullfile(out_dir,['phz_Esens_it',num2str(Itr-1),tag,'.fits']));
fitswrite_tb(mp,tb,Im_tb.Iinco,fullfile(out_dir,['normI_inco_it',num2str(Itr-1),tag,'.fits']));

if(~strcmpi(mp.estimator,'perfect'))
    ev = Im_tb.ev;
    save(fullfile(out_dir,['probing_data_',num2str(Itr-1),tag,'.mat']),'ev');
end
if strcmpi(mp.estimator,'ekf_maintenance') && any(mp.est.itr_ol==ev.Itr) == true
    img = mean(Im_tb.ev.normI_OL_sbp,3);
    fitswrite_tb(mp,tb,img,fullfile(out_dir,['normI_OL_it',num2str(Itr-1),tag,'.fits']));
end
% Update the diary 
diary off; diary(mp.diaryfile)

end %--END OF FUNCTION

function fitswrite_tb(mp, tb, obj, filename)
if mp.flagSim
    fitswrite(obj,filename);
else
    sciCam_fitswrite(tb,obj,filename);
end
end