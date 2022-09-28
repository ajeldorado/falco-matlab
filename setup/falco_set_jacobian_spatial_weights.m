% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Assign weights to the control Jacobian by spatial location.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_set_jacobian_spatial_weights(mp)

    [XISLAMD,ETASLAMD] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
    RHOScompact = sqrt(XISLAMD.^2+ETASLAMD.^2);

    %--Spatial weighting matrix
    mp.Wspatial = mp.Fend.corr.mask;
    if ~isempty(mp.WspatialDef)
        for kk=1:size(mp.WspatialDef,1)
            Wannulus = 1+(sqrt(mp.WspatialDef(kk,3))-1)*((RHOScompact>=mp.WspatialDef(kk,1)) & (RHOScompact<mp.WspatialDef(kk,2)));
            mp.Wspatial = mp.Wspatial.*Wannulus;
        end
    end

    %--Spatial weighting vector (for each star)
    Npix = sum(sum(mp.Fend.corr.maskBool));
%     mp.WspatialVec = zeros(Npix, mp.jac.Nmode);
%     for iMode = 1:mp.jac.Nmode
%         iStar = mp.jac.star_inds(iMode);
%         mp.WspatialVec(:, iMode) = mp.jac.star.weights(iStar) * mp.Wspatial(mp.Fend.corr.maskBool); 
%     end
    mp.WspatialVec = zeros(Npix, mp.compact.star.count);
    for iStar = 1:mp.compact.star.count
        mp.WspatialVec(:, iStar) = mp.jac.star.weights(iStar) * mp.Wspatial(mp.Fend.corr.maskBool); 
    end
    if strcmpi(mp.estimator, 'iefc')
        mp.WspatialVec = repmat(mp.WspatialVec, [size(mp.iefc.probeCube, 3), 1]);
    end
    
%     mp.WspatialVec = mp.Wspatial(mp.Fend.corr.maskBool); 
    if(mp.flagFiber);  mp.WspatialVec = ones(mp.Fend.Nfiber, 1);  end
    if(mp.flagFiber && mp.flagLenslet);  mp.WspatialVec = ones(mp.Fend.Nlens, 1);  end

end
