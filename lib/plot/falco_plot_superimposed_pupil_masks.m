% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Plot the pupil masks on top of each other to make sure they are oriented
% correctly.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% None

function falco_plot_superimposed_pupil_masks(mp)

    %--Plot superposed apodizer and telescope pupil
    if mp.flagApod && mp.flagPlot
        figure(600)
        imagesc(mp.P3.compact.mask - mp.P1.compact.mask,[-1, 1]);
        title('Apodizer - Entrance Pupil', 'Fontsize', 16)
        axis xy equal tight; colorbar;
        set(gca, 'Fontsize', 20);
        set(gcf, 'Color', 'w');
        drawnow
    end
    
    %--Plot superposed Lyot stop and telescope pupil
    switch upper(mp.coro)
        case{'HLC', 'LC', 'APLC', 'VC', 'AVC', 'VORTEX'}
            
            if mp.flagPlot
                P4mask = pad_crop(mp.P4.compact.mask,mp.P1.compact.Narr);
                P4mask = propcustom_relay(P4mask, mp.Nrelay1to2+mp.Nrelay2to3+mp.Nrelay3to4);
                P1andP4 = abs(mp.P1.compact.mask) + P4mask;
                
                figure(301)
                imagesc(P1andP4);
                title('Superimposed Pupil and Lyot Stop', 'Fontsize', 16)
                axis xy equal tight; colorbar; 
                set(gca, 'Fontsize', 20);
                set(gcf, 'Color', 'w');
                drawnow;

                if mp.flagApod
                    P1andP3 = abs(mp.P1.compact.mask) + ...
                        pad_crop(mp.P3.compact.mask,length(mp.P1.compact.mask));
                    
                    figure(302)
                    imagesc(P1andP3);
                    title('Superimposed Pupil and Apodizer', 'Fontsize', 16);
                    axis xy equal tight; colorbar; 
                    set(gca, 'Fontsize', 20);
                    set(gcf, 'Color', 'w');
                    drawnow;
                end
            end
    end
    
end
