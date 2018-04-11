% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function dm = falco_gen_dm_surf(dm)
% 
%--Function to generate a deformable mirror (DM) surface using PROPER.
%
%--INPUTS
% dm: structure of DM parameters
%
%--OUTPUT
% DMsurf: DM surface in meters
% 
%--VERSION HISTORY
% -Modified on 2018-02-06 by A.J. Riggs to have dx and N be inputs as well 
% in order to avoid having to define a few variables (such as the DM 
% commands) twice for the compact and full models.
% -Created falco_gen_dm_poke_cube_PROPER.m on 2017-11-17 by A.J. Riggs. 
%
function DMsurf = falco_gen_dm_surf(dm,dx,N)

%--Set the order of operations
orderOfOps = 'XYZ';
if(isfield(dm,'flagZYX')); 
    if(dm.flagZYX)
        orderOfOps = 'ZYX'; 
    end
end

pupil_ratio = 1;
wl_dummy = 1e-6; %--dummy value needed to initialize wavelength in PROPER (meters)

bm  = prop_begin(N*dx, wl_dummy, N, pupil_ratio);
[~,DMsurf] = prop_dm(bm, dm.VtoH.*dm.V, dm.xc, dm.yc, dm.dm_spacing,'XTILT',dm.xtilt,'YTILT',dm.ytilt,'ZTILT',dm.zrot,orderOfOps);
%     figure(1); imagesc(DMsurf); axis xy equal tight; colorbar; 

end %--END OF FUNCTION

% function DMsurf = falco_gen_dm_surf(dm)
% 
% %--Set the order of operations
% orderOfOps = 'XYZ';
% if(isfield(dm,'flagZYX')); 
%     if(dm.flagZYX)
%         orderOfOps = 'ZYX'; 
%     end
% end
% 
% pupil_ratio = 1;
% wl_dummy = 1e-6; %--dummy value needed to initialize wavelength in PROPER (meters)
% 
% bm  = prop_begin(dm.NdmPad*dm.dx, wl_dummy, dm.NdmPad, pupil_ratio);
% [~,DMsurf] = prop_dm(bm, dm.VtoH.*dm.V, dm.xc, dm.yc, dm.dm_spacing,'XTILT',dm.xtilt,'YTILT',dm.ytilt,'ZTILT',dm.zrot,orderOfOps);
% %     figure(1); imagesc(DMsurf); axis xy equal tight; colorbar; 
% 
% end %--END OF FUNCTION







