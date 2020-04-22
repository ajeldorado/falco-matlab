% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% REVISION HISTORY:
% --------------
% Created by G. Ruane on 2020-04-21
% ---------------

function mp = falco_configure_WFS(mp)

    if(mp.flagWFS)
        
        mp = falco_gen_FPM_ZWFS(mp);
        mp = falco_get_FPM_ZWFS_coordinates(mp);
        mp = falco_set_WFS_spectral_properties(mp);
        
        
        % Intialize E-field seen by WFS (will be different to coronagraph 
        % if using WFS out of band). This is updated in model_ZWFS.
        mp.wfs.E = ones(mp.P1.full.Narr,mp.P1.full.Narr,mp.wfs.Nwpsbp,mp.wfs.Nsbp);
        
    end
        
end %--END OF FUNCTION