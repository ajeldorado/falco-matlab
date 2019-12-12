% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%
% REVISION HISTORY:
% --------------
% Created by A.J. Riggs on 2018-10-01 by extracting material from
% falco_init_ws.m.
% ---------------

function [mp] = falco_gen_FPM_LC(mp)

    if(mp.compact.flagGenFPM || mp.full.flagGenFPM)
        %--Make or read in focal plane mask (FPM) amplitude for the full model
        FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
        FPMgenInputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
        FPMgenInputs.FPMampFac = mp.FPMampFac; % amplitude transmission of inner FPM spot
        FPMgenInputs.centering = mp.centering;
    end
        
    if(mp.full.flagGenFPM)
        FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
        mp.F3.full.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
    end
    mp.F3.full.Nxi = size(mp.F3.full.mask.amp,2);
    mp.F3.full.Neta= size(mp.F3.full.mask.amp,1);   

    
    if(mp.compact.flagGenFPM)
        %--Number of points across the FPM in the compact model
        if(isinf(mp.F3.Rout))
            switch mp.centering
            case 'pixel'
                mp.F3.compact.Nxi = ceil_even((2*(mp.F3.Rin*mp.F3.compact.res + 1/2)));
            case 'interpixel'
                mp.F3.compact.Nxi = ceil_even((2*mp.F3.Rin*mp.F3.compact.res));
            end
        else
            switch mp.centering
                case 'pixel'
                    mp.F3.compact.Nxi = ceil_even((2*(mp.F3.Rout*mp.F3.compact.res + 1/2)));
                case 'interpixel'
                    mp.F3.compact.Nxi = ceil_even((2*mp.F3.Rout*mp.F3.compact.res));
            end
        end
        mp.F3.compact.Neta = mp.F3.compact.Nxi;
        
        %--Make or read in focal plane mask (FPM) amplitude for the compact model
        FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
        mp.F3.compact.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
    else
        mp.F3.compact.Nxi = size(mp.F3.compact.mask.amp,2);
        mp.F3.compact.Neta= size(mp.F3.compact.mask.amp,1); 
    end

end %--END OF FUNCTION