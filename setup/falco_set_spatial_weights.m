% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Spatial weighting of pixel intensity in the dark hole

function mp = falco_set_spatial_weights(mp)

    [XISLAMD,ETASLAMD] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
    RHOScompact = sqrt(XISLAMD.^2+ETASLAMD.^2);

    %--Spatial weighting matrix
    mp.Wspatial = mp.Fend.corr.mask;
    if(isempty(mp.WspatialDef)==false)
        for kk=1:size(mp.WspatialDef,1)
            Wannulus = 1+(sqrt(mp.WspatialDef(kk,3))-1)*((RHOScompact>=mp.WspatialDef(kk,1)) & (RHOScompact<mp.WspatialDef(kk,2)));
            mp.Wspatial = mp.Wspatial.*Wannulus;
        end
    end

    %--Spatial weighting vector
    mp.WspatialVec = mp.Wspatial(mp.Fend.corr.maskBool); 
    if(mp.flagFiber && mp.flagLenslet);  mp.WspatialVec = ones(mp.Fend.Nlens,1);  end

end