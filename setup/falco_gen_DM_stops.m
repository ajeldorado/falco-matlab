% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate circular aperture stops at DMs 1 and 2.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_gen_DM_stops(mp)

    if mp.flagDM1stop
        mp.dm1.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm1.Dstop,mp.centering);
        mp.dm1.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm1.Dstop,mp.centering);
    end
    
    if mp.flagDM2stop
        mp.dm2.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm2.Dstop,mp.centering);
        mp.dm2.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm2.Dstop,mp.centering);
    end

end
