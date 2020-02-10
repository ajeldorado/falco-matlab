% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% FPM coordinates, [meters] and [dimensionless]

function mp = falco_get_FPM_coordinates(mp)


switch upper(mp.coro)
    case{'VORTEX','VC','AVC'}   %--Nothing needed to run the vortex model
    case 'SPHLC' %--Moved to separate function
    otherwise
        
        switch mp.layout
            case{'wfirst_phaseb_simple','wfirst_phaseb_proper'}
            otherwise
                %--FPM (at F3) Resolution [meters]
                mp.F3.full.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.res;
                mp.F3.full.deta = mp.F3.full.dxi;
        end
        mp.F3.compact.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.res;
        mp.F3.compact.deta = mp.F3.compact.dxi;
        
        %--Coordinates in FPM plane in the compact model [meters]
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.compact.Nxi,2)==1)
            mp.F3.compact.xis  = (-(mp.F3.compact.Nxi-1)/2:(mp.F3.compact.Nxi-1)/2)*mp.F3.compact.dxi;
            mp.F3.compact.etas = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2).'*mp.F3.compact.deta;
        else
            mp.F3.compact.xis  = (-mp.F3.compact.Nxi/2: (mp.F3.compact.Nxi/2-1))*mp.F3.compact.dxi;
            mp.F3.compact.etas = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1)).'*mp.F3.compact.deta;
        end

        switch mp.layout
            case{'wfirst_phaseb_simple','wfirst_phaseb_proper'}
            otherwise
                %--Coordinates (dimensionless [DL]) for the FPMs in the full model
                if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.full.Nxi,2)==1)
                    mp.F3.full.xisDL  = (-(mp.F3.full.Nxi-1)/2:(mp.F3.full.Nxi-1)/2)/mp.F3.full.res;
                    mp.F3.full.etasDL = (-(mp.F3.full.Neta-1)/2:(mp.F3.full.Neta-1)/2)/mp.F3.full.res;
                else
                    mp.F3.full.xisDL  = (-mp.F3.full.Nxi/2:(mp.F3.full.Nxi/2-1))/mp.F3.full.res;
                    mp.F3.full.etasDL = (-mp.F3.full.Neta/2:(mp.F3.full.Neta/2-1))/mp.F3.full.res;
                end
        end
        
        %--Coordinates (dimensionless [DL]) for the FPMs in the compact model
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.compact.Nxi,2)==1)
            mp.F3.compact.xisDL  = (-(mp.F3.compact.Nxi-1)/2:(mp.F3.compact.Nxi-1)/2)/mp.F3.compact.res;
            mp.F3.compact.etasDL = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2)/mp.F3.compact.res;
        else
            mp.F3.compact.xisDL  = (-mp.F3.compact.Nxi/2:(mp.F3.compact.Nxi/2-1))/mp.F3.compact.res;
            mp.F3.compact.etasDL = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1))/mp.F3.compact.res;
        end
end


end %--END OF FUNCTION