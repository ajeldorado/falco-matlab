% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Plot relevant data for Borde-Traub probing.

function falco_plot_borde_traub_probes(mp, ev, dDMVplus, iSubband)

if mp.flagPlot
    
    if mp.est.probe.whichDM == 1
        VtoH = mp.dm1.VtoH;
    elseif mp.est.probe.whichDM == 2
        VtoH = mp.dm2.VtoH;
    end
    
    Npairs = mp.est.probe.Npairs;

    subplot = @(m,n,p) subtightplot(m,n,p,[0.025 0.025],[0.1 0.1],[0.1 0.1]); %--Define a trimmed-down plot

    handles.master = figure(90+ev.iStar);
    set(handles.master,'units', 'inches', 'Position', [0 0 12 8])
    set(handles.master,'Color','w')

    %--Plot the DM shapes for each probe
    for iProbe=1:Npairs
        subplot(2,Npairs,iProbe); % Save the handle of the subplot
        imagesc(1e9*dDMVplus(:,:,iProbe).*VtoH); axis xy equal tight; axis off;
        colorbar;
        colormap(gca,parula);
        title(sprintf('DM Probe %d (nm)',iProbe));
    end
    hold on;
    

    %--Plot the raw images for each + probe
    IcubeNonNeg = ev.imageArray(:, :, :, iSubband);
    IcubeNonNeg(IcubeNonNeg<0) = 0;
    for iProbe=1:Npairs
        subplot(2,Npairs,iProbe + Npairs*1); % Save the handle of the subplot
        imagesc(log10(IcubeNonNeg(:,:,iProbe*2))); axis xy equal tight; axis off;
        colorbar;
        colormap(gca,gray);
        title(sprintf('+ Probe Image'));
    end
    
    drawnow;
    hold off

end