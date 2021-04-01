% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

function [mp] = falco_gen_FPM_LC(mp)

    if(mp.compact.flagGenFPM || mp.full.flagGenFPM)
        %--Make or read in focal plane mask (FPM) amplitude for the full model
        FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
        FPMgenInputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
        FPMgenInputs.centering = mp.centering;
        if isfield(mp, 'FPMampFac'); FPMgenInputs.FPMampFac = mp.FPMampFac; end % amplitude transmission of inner FPM spot
    end
        
    
    if(mp.full.flagGenFPM)
        FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
        if(isfield(mp.F3.full,'mask'))
            mp.F3.full = rmfield(mp.F3.full,'mask');
        end
        mp.F3.full.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
    end
 

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
        if(isfield(mp.F3.compact,'mask'))
            mp.F3.compact = rmfield(mp.F3.compact,'mask');
        end
        mp.F3.compact.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
    else
        %
    end

end %--END OF FUNCTION