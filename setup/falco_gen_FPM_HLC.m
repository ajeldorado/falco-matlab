% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% REVISION HISTORY:
% --------------
% Created by A.J. Riggs on 2018-10-01 by extracting material from
% falco_init_ws.m.
% ---------------

function [mp] = falco_gen_FPM_HLC(mp)

    %% Coordinates

    %--Number of points across the FPM in the full model
    mp.F3.full.Nxi  = mp.dm9.NdmPad;
    mp.F3.full.Neta = mp.F3.full.Nxi;
    %--Number of points across the FPM in the compact model
    mp.F3.compact.Nxi  = mp.dm9.compact.NdmPad;
    mp.F3.compact.Neta = mp.F3.compact.Nxi;

    %% Starting DM8 and DM9 surfaces
    
        if(mp.flagPlot)
            DM8surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm8,'full');
            figure(8); imagesc(DM8surf); axis xy equal tight; colorbar; title('FPM Metal Thickness','Fontsize',20); set(gca,'Fontsize',20); drawnow; 
        end
        if(mp.flagPlot)
            DM9surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm9,'full');
            figure(9); imagesc(DM9surf); axis xy equal tight; colorbar; title('FPM Dielectric Thickness','Fontsize',20); set(gca,'Fontsize',20); drawnow; 
        end

%         for ii = 1:mp.dm9.NactTotal
%             mp.dm9.V = 0*mp.dm9.V;
%             mp.dm9.V(ii) = 1;
%             DM9surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm9, 'full');
%             figure(201); imagesc(DM9surf); axis xy equal tight; colorbar; title(sprintf('%04d', ii), 'Fontsize', 20); drawnow;
%             pause(1/25);
%         end

end %--END OF FUNCTION