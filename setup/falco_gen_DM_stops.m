% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% 
%  (moved here because the commands mp.dm2.compact = mp.dm2; and mp.dm1.compact = mp.dm1; otherwise would overwrite the compact model masks)

function [mp] = falco_gen_DM_stops(mp)

    %%--DM Aperture Masks 
    if(mp.flagDM1stop)
        mp.dm1.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm1.Dstop,mp.centering);
        mp.dm1.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm1.Dstop,mp.centering);
    end
    
    if(mp.flagDM2stop)
        mp.dm2.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm2.Dstop,mp.centering);
        mp.dm2.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm2.Dstop,mp.centering);
    end

end %--END OF FUNCTION