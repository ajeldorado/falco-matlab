% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% FPM coordinates, [meters] and [dimensionless]

function mp = falco_get_FPM_ZWFS_coordinates(mp)

    lam0 = mp.wfs.lambda0;

    mp.wfs.mask.dxi = (mp.fl*lam0/mp.P2.D)/mp.wfs.mask.res;
    mp.wfs.mask.deta = mp.wfs.mask.dxi;

    %--Coordinates in FPM plane in the compact model [meters]
    if(strcmpi(mp.centering,'interpixel') || mod(mp.wfs.mask.Nxi,2)==1)
        mp.wfs.mask.xis  = (-(mp.wfs.mask.Nxi-1)/2:(mp.wfs.mask.Nxi-1)/2)*mp.wfs.mask.dxi;
        mp.wfs.mask.etas = (-(mp.wfs.mask.Neta-1)/2:(mp.wfs.mask.Neta-1)/2).'*mp.wfs.mask.deta;
    else
        mp.wfs.mask.xis  = (-mp.wfs.mask.Nxi/2: (mp.wfs.mask.Nxi/2-1))*mp.wfs.mask.dxi;
        mp.wfs.mask.etas = (-mp.wfs.mask.Neta/2:(mp.wfs.mask.Neta/2-1)).'*mp.wfs.mask.deta;
    end

    %--Coordinates (dimensionless [DL]) for the FPMs in the compact model
    if(strcmpi(mp.centering,'interpixel') || mod(mp.wfs.mask.Nxi,2)==1)
        mp.wfs.mask.xisDL  = (-(mp.wfs.mask.Nxi-1)/2:(mp.wfs.mask.Nxi-1)/2)/mp.wfs.mask.res;
        mp.wfs.mask.etasDL = (-(mp.wfs.mask.Neta-1)/2:(mp.wfs.mask.Neta-1)/2)/mp.wfs.mask.res;
    else
        mp.wfs.mask.xisDL  = (-mp.wfs.mask.Nxi/2:(mp.wfs.mask.Nxi/2-1))/mp.wfs.mask.res;
        mp.wfs.mask.etasDL = (-mp.wfs.mask.Neta/2:(mp.wfs.mask.Neta/2-1))/mp.wfs.mask.res;
    end


end %--END OF FUNCTION