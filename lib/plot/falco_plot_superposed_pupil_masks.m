% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Plot the pupil and Lyot stop on top of each other to make sure they are aligned correctly
% 
% Only for coronagraphs using Babinet's principle, for which the input pupil and Lyot plane have 
% the same resolution.
function falco_plot_superposed_pupil_masks(mp)

    %--Plot superposed apodizer and telescope pupil
    if(mp.flagApod && mp.flagPlot)
        figure(600); imagesc(mp.P3.compact.mask - mp.P1.compact.mask,[-1 1]); axis xy equal tight; colorbar; drawnow;
    end
    
    %--Plot superposed Lyot stop and telescope pupil
    switch upper(mp.coro)
        case{'FOHLC','HLC','LC','APLC','VC','AVC'}
            if(mp.flagPlot)
                P4mask = padOrCropEven(mp.P4.compact.mask,mp.P1.compact.Narr);
                P4mask = propcustom_relay(P4mask, mp.Nrelay1to2+mp.Nrelay2to3+mp.Nrelay3to4);
%                 P4mask = rot90(P4mask,2);
%                 if(strcmpi(mp.centering,'pixel'))
%                    P4mask = circshift(P4mask,[1 1]); 
%                 end
                P1andP4 = mp.P1.compact.mask + P4mask;
                figure(301); imagesc(P1andP4); axis xy equal tight; colorbar; 
                set(gca,'Fontsize',20); title('Pupil and LS Superimposed','Fontsize',16');

                if(mp.flagApod)
                    P1andP3 = mp.P1.compact.mask + ...
                        padOrCropEven(mp.P3.compact.mask,length(mp.P1.compact.mask));
                    figure(302); imagesc(P1andP3); axis xy equal tight; colorbar; 
                    set(gca,'Fontsize',20); title('Pupil and Apod Superimposed','Fontsize',16');
                end
            end
    end
    
end