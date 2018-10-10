% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Spatial weighting of pixel intensity. 
% NOTE: For real instruments and testbeds, only the compact model should be 
% used. The full model spatial weighting is included too if in simulation 
% the full model has a different detector resolution than the compact model.
%
% REVISION HISTORY:
% --------------
% Created by A.J. Riggs on 2018-10-01 by extracting material from falco_init_ws.m.
% ---------------

function mp = falco_config_spatial_weights(mp)



[XISLAMD,ETASLAMD] = meshgrid(mp.F4.xisDL, mp.F4.etasDL);
RHOScompact = sqrt(XISLAMD.^2+ETASLAMD.^2);
% [XISLAMD,ETASLAMD] = meshgrid(mp.F4.xisDL, mp.F4.etasDL);
% RHOSfull = sqrt(XISLAMD.^2+ETASLAMD.^2);

mp.Wspatial = mp.F4.corr.mask;
% mp.WspatialFull = mp.F4.corr.mask;
if( isfield(mp,'WspatialDef') ) %--If there are spatial weights for the Jacobian
    if(isempty(mp.WspatialDef)==false)
    
        for kk=1:size(mp.WspatialDef,1)
            Wannulus = 1+(sqrt(mp.WspatialDef(kk,3))-1)*( (RHOScompact>=mp.WspatialDef(kk,1)) & (RHOScompact<mp.WspatialDef(kk,2)) );
            %WannulusFull = 1+(sqrt(mp.WspatialDef(kk,3))-1)*( (RHOSfull>=mp.WspatialDef(kk,1)) & (RHOSfull<mp.WspatialDef(kk,2)) );
            %if(mp.flagPlot); figure(70); imagesc(mp.F4.xisDL,mp.F4.etasDL,WannulusCompact); axis xy equal tight; colorbar; pause(1); end
            
            mp.Wspatial = mp.Wspatial.*Wannulus;
            %mp.WspatialFull = mp.WspatialFull.*WannulusFull;
            %if(mp.flagPlot); figure(71); imagesc(mp.F4.xisDL,mp.F4.etasDL,mp.Wspatial); axis xy equal tight; colorbar; drawnow; end
        end
    else
        %--nothing
    end
    
else
 %--nothing
end
