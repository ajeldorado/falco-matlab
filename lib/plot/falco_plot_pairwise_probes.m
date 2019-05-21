% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to plot relevant data for pairwise probing.
%
% REVISION HISTORY:
% -Created on 2018-04-22 by A.J. Riggs.

function falco_plot_pairwise_probes(mp,ev,dDMVplus,ampSq2Dcube)

if(mp.flagPlot)
    Npairs = mp.est.probe.Npairs;

    subplot = @(m,n,p) subtightplot(m,n,p,[0.025 0.025],[0.1 0.1],[0.1 0.1]); %--Define a trimmed-down plot

    handles.master = figure(99);
    set(handles.master,'units', 'inches', 'Position', [0 0 12 8])
    set(handles.master,'Color','w')

    %--Plot the DM shapes for each probe
    for iProbe=1:Npairs
        subplot(4,Npairs,iProbe); % Save the handle of the subplot
        imagesc(1e9*dDMVplus(:,:,iProbe).*mp.dm1.VtoH); axis xy equal tight; axis off;
        colorbar;
        colormap(gca,parula);
        title(sprintf('Probe %d (nm)',iProbe));
    end
    hold on;
    

    %--Plot the raw images for each + probe
    IcubeNonNeg = ev.Icube;
    IcubeNonNeg(IcubeNonNeg<0) = 0;
    for iProbe=1:Npairs
        subplot(4,Npairs,iProbe + Npairs*1); % Save the handle of the subplot
        imagesc(log10(IcubeNonNeg(:,:,iProbe*2))); axis xy equal tight; axis off;
        colorbar;
        colormap(gca,gray);
        title(sprintf('+ Probe Image'));
    end
    
    %--Plot the raw images for each - probe
    IcubeNonNeg = ev.Icube;
    IcubeNonNeg(IcubeNonNeg<0) = 0;
    for iProbe=1:Npairs
        subplot(4,Npairs,iProbe + Npairs*2); % Save the handle of the subplot
        imagesc(log10(IcubeNonNeg(:,:,iProbe*2+1))); axis xy equal tight; axis off;
        colorbar;
        colormap(gca,gray);
        title(sprintf('- Probe Image'));
    end
    
    %--Plot the squared amplitude of the delta amplitude from each probe.
    for iProbe=1:Npairs
        subplot(4,Npairs,iProbe + Npairs*3); % Save the handle of the subplot
        ampSq2D = ampSq2Dcube(:,:,iProbe);
        imagesc(log10(ampSq2D),log10([max(ampSq2D(:))*1e-1,max(ampSq2D(:))])); axis xy equal tight; axis off;
        colorbar;
        colormap(gca,gray);
        title('Probe Intensity, |dP|^2');
    end
    drawnow;
    hold off

end