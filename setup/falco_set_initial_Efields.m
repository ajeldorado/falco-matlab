% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Initial Electric Fields for Star and Exoplanet
% (moved here because the commands mp.dm2.compact = mp.dm2; and mp.dm1.compact = mp.dm1; otherwise would overwrite the compact model masks)

function [mp] = falco_set_initial_Efields(mp)

    if(isfield(mp.P1.full,'E')==false)
        mp.P1.full.E  = ones(mp.P1.full.Narr,mp.P1.full.Narr,mp.Nwpsbp,mp.Nsbp);
    end
    
    if(isfield(mp,'Eplanet')==false)
        mp.Eplanet = mp.P1.full.E; % NOTE: Phase ramp added later in propagation model
    end
    
    if(isfield(mp.P1.compact,'E')==false)
        mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp);
    end

end %--END OF FUNCTION