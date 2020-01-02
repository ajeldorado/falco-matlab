% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%

function mp = falco_set_dm_surface_padding(mp)

    %--Compact Model: Set nominal DM plane array sizes as a power of 2 for angular spectrum propagation with FFTs
    if(any(mp.dm_ind==1) && any(mp.dm_ind==2))
        NdmPad = 2.^ceil(1 + log2(max([mp.dm1.compact.NdmPad,mp.dm2.compact.NdmPad]))); 
    elseif(any(mp.dm_ind==1))
        NdmPad = 2.^ceil(1 + log2(mp.dm1.compact.NdmPad));
    elseif(any(mp.dm_ind==2))
        NdmPad = 2.^ceil(1 + log2(mp.dm2.compact.NdmPad));
    else
        NdmPad = 2*mp.P1.compact.Nbeam;
    end
    while((NdmPad < min(mp.sbp_centers)*abs(mp.d_dm1_dm2)/mp.P2.full.dx^2) || (NdmPad < min(mp.sbp_centers)*abs(mp.d_P2_dm1)/mp.P2.compact.dx^2)) %--Double the zero-padding until the angular spectrum sampling requirement is not violated
        NdmPad = 2*NdmPad; 
    end
    mp.compact.NdmPad = NdmPad;

    %--Full Model: Set nominal DM plane array sizes as a power of 2 for angular spectrum propagation with FFTs
    if(any(mp.dm_ind==1) && any(mp.dm_ind==2))
        NdmPad = 2.^ceil(1 + log2(max([mp.dm1.NdmPad,mp.dm2.NdmPad]))); 
    elseif(any(mp.dm_ind==1))
        NdmPad = 2.^ceil(1 + log2(mp.dm1.NdmPad));
    elseif(any(mp.dm_ind==2))
        NdmPad = 2.^ceil(1 + log2(mp.dm2.NdmPad));
    else
        NdmPad = 2*mp.P1.full.Nbeam;    
    end
    while((NdmPad < min(mp.full.lambdas)*abs(mp.d_dm1_dm2)/mp.P2.full.dx^2) || (NdmPad < min(mp.full.lambdas)*abs(mp.d_P2_dm1)/mp.P2.full.dx^2)) %--Double the zero-padding until the angular spectrum sampling requirement is not violated
        NdmPad = 2*NdmPad; 
    end
    mp.full.NdmPad = NdmPad;

end