% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%
% REVISION HISTORY:
% --------------
% Created by A.J. Riggs on 2018-10-01 by extracting material from
% falco_init_ws.m.
% ---------------

function [mp] = falco_config_gen_FPM_HLC(mp)


    %% Coordinates

    %--Number of points across the FPM in the full model
    mp.F3.full.Nxi  = mp.dm9.NdmPad;%ceil_even( 2*(mp.F3.Rin+buffer)*mp.F3.full.res);
    mp.F3.full.Neta = mp.F3.full.Nxi;
    %--Number of points across the FPM in the compact model
    mp.F3.compact.Nxi  = mp.dm9.compact.NdmPad; %ceil_even( 2*(mp.F3.Rin+buffer)*mp.F3.compact.res);
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


end %--END OF FUNCTION
