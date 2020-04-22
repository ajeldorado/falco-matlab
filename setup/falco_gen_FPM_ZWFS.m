% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% REVISION HISTORY:
% --------------
% Created by A.J. Riggs on 2018-10-01 by extracting material from
% falco_init_ws.m.
% Modifed from falco_config_gen_FPM_LC by G. Ruane on 2018-11-20.
% Modifed from falco_gen_FPM_Roddier by G. Ruane on 2020-04-21.
% ---------------

function mp = falco_gen_FPM_ZWFS(mp)

        %--Make focal plane mask (FPM) amplitude
        FPMgenInputs.pixresFPM = mp.wfs.mask.res; %--pixels per lambda_c/D
        FPMgenInputs.rhoInner = mp.wfs.mask.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
        FPMgenInputs.rhoOuter = mp.wfs.mask.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
        FPMgenInputs.FPMampFac = mp.wfs.mask.FPMampFac; % amplitude transmission of inner FPM spot
        FPMgenInputs.centering = mp.centering;
        mp.wfs.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
        
        %--Number of points across the FPM
        if(isinf(mp.wfs.mask.Rout))
            switch mp.centering
            case 'pixel'
                mp.wfs.mask.Nxi = ceil_even((2*(mp.wfs.mask.Rin*mp.wfs.mask.res + 1/2)));
            case 'interpixel'
                mp.wfs.mask.Nxi = ceil_even((2*mp.wfs.mask.Rin*mp.wfs.mask.res));
            end
        else
            switch mp.centering
                case 'pixel'
                    mp.wfs.mask.Nxi = ceil_even((2*(mp.wfs.mask.Rout*mp.wfs.mask.res + 1/2)));
                case 'interpixel'
                    mp.wfs.mask.Nxi = ceil_even((2*mp.wfs.mask.Rout*mp.wfs.mask.res));
            end
        end
        mp.wfs.mask.Neta = mp.wfs.mask.Nxi;
        
        %-- Make a mask for the phase pattern support
        FPMgenInputsPhz = FPMgenInputs;
        FPMgenInputsPhz.FPMampFac = 0;
        mp.wfs.mask.phzSupport = 1-falco_gen_annular_FPM(FPMgenInputsPhz);
        
        % need to define the index of refraction in the case of a
        % transmissive mask 
        if(strcmpi(mp.wfs.mask.type,'transmissive'))
            if(strcmp(mp.wfs.mask.material,'FS'))% Fused silica 
                nsquared = @(lam_um) 1 + (0.6961663*lam_um^2)/(lam_um^2-0.0684043^2)+ ...
                                      (0.4079426*lam_um^2)/(lam_um^2-0.1162414^2)+ ...
                                      (0.8974794*lam_um^2)/(lam_um^2-9.896161^2);          
                n = @(lam) real(sqrt(nsquared(lam*1e6)));
                mp.wfs.mask.n = n;
            else
                error('Material not defined for Zernike WFS mask.');
            end
        end

end %--END OF FUNCTION