% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function handles = falco_plot_progress_gpct(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf)

subplot = @(m,n,p) subtightplot(m,n,p,[0.025 0.025],[0.1 0.1],[0.1 0.1]);

Icbmin = -8;
Icbmax = -4;
Im = Im_tb.Im;
if(Itr>1)
    Imod = Inorm.mod(Itr-1);
else
    Imod = NaN;
end
Im(Im<0) = 0; %--Prevent the log10(Im) plot from getting complex values.

if(mp.flagPlot)

    if(Itr>1)
        delete(handles.tb1)
        delete(handles.tb2)
        delete(handles.tb3)
        delete(handles.tb4)
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

    
    subplot(2,3,1); % Save the handle of the subplot
    axis off
    handles.tb1 = text(0.1,0.8,sprintf('%s: Iteration %d',mp.coro,Itr-1));
    handles.tb2 = text(0.1,0.7,sprintf('%.1f%% BW @ %dnm',(100*mp.fracBW),round(mp.lambda0*1e9)));
    handles.tb3 = text(0.1,0.6,sprintf('I_{norm} = %.2e',Inorm.total(Itr)));
    handles.tb4 = text(0.1,0.5,sprintf('I_{mod,prev} = %.2e',Imod));
    switch lower(mp.thput_metric)
        case{'hmi'} %--Absolute energy within half-max isophote(s)
            handles.tb5 = text(0.1,0.4,sprintf('T_{half-max} =   %.2f%%',100*mp.thput_vec(Itr)));
        case{'ee','e.e.'} %--Absolute energy encircled within a given radius
            handles.tb5 = text(0.1,0.4,sprintf('T_{E.E.} =   %.2f%%',100*mp.thput_vec(Itr)));
    end

    subplot(2,3,2); % Save the handle of the subplot
    imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im),[Icbmin Icbmax]); 
    axis xy equal tight; 
    colorbar; 
    colormap(gca,parula);
%     xlabel('\lambda_0/D'); 
%     ylabel('\lambda_0/D');
%     ylabel(ch_psf,'log(NI)');
    title('Stellar PSF');
%     title(sprintf('PSF, NI=%.2e',contrast_bandavg(Itr)),'Fontsize',20,'Fontweight','Bold');
%       title(sprintf('PSF at Iter=%03d,   NI=%.2e',Itr,contrast_bandavg(Itr)),'Fontsize',20,'Fontweight','Bold');


	subplot(2,3,3); % Save the handle of the subplot
    imagesc(1e9*DM1surf);  axis xy equal tight; axis off;
    colorbar;
    colormap(gca,gray);
    title('DM1 Surface (nm)');

    subplot(2,3,4);
    semilogy(0:length(Inorm.total)-1,Inorm.total,'-o');hold on;
    if(Itr>1)
        semilogy(0:Itr-2,Inorm.mod,'-o');
        semilogy(0:Itr-2,Inorm.unmod,'--o');
    end
    hold off;
    xlim([0 length(Inorm.total)])
    xlabel('Iteration')
%     ylabel('Norm. I');
    legend('Total','Modulated','Unmodulated');
	title('Normalized Intensity')
    grid on;axis square;
% 	hcbdummy = colorbar;set(hcbdummy,'visible','off');
    
	subplot(2,3,5); % Save the handle of the subplot
    imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(abs(Im_tb.E).^2),[Icbmin Icbmax]); 
    axis xy equal tight;
    colorbar;
    colormap(gca,parula)
%     xlabel('\lambda_0/D'); 
%     ylabel('\lambda_0/D');
    title('Modulated (previous)');
    
	subplot(2,3,6); % Save the handle of the subplot
    imagesc(mp.Fend.xisDL,mp.Fend.etasDL,angle(Im_tb.E),[-pi pi]); 
    axis xy equal tight; 
    colorbar; 
    colormap(gca,hsv);
%     xlabel('\lambda_0/D'); 
%     ylabel('\lambda_0/D');
    title('Phase (previous)');
    
   drawnow;
   
end

out_dir = [mp.bench.info.OUT_DATA_DIR,mp.runLabel,'/'];
% Directory to save dat
if(~exist(out_dir, 'dir'))
    mkdir(out_dir);
end

hcst_andor_fitswrite(mp.bench,Im,[out_dir,'normI_it',num2str(Itr-1),'.fits'],false);
hcst_andor_fitswrite(mp.bench,mp.dm1.V,[out_dir,'dmV_it',num2str(Itr-1),'.fits'],false);
hcst_andor_fitswrite(mp.bench,DM1surf,[out_dir,'dmmodel_it',num2str(Itr-1),'.fits'],false);

if(Itr>1)
    hcst_andor_fitswrite(mp.bench,abs(Im_tb.E).^2,[out_dir,'normI_Esens_it',num2str(Itr-2),'.fits'],false);
    hcst_andor_fitswrite(mp.bench,angle(Im_tb.E),[out_dir,'phz_Esens_it',num2str(Itr-2),'.fits'],false);
end

end %--END OF FUNCTION
