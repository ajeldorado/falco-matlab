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

function [mp] = falco_config_gen_FPM_LC(mp)

        %--Make or read in focal plane mask (FPM) amplitude for the full model
        FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
        FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
        FPMgenInputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
        FPMgenInputs.FPMampFac = mp.FPMampFac; % amplitude transmission of inner FPM spot
        FPMgenInputs.centering = mp.centering;
        mp.F3.full.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
        % figure(204); imagesc(mp.F3.full.mask.amp); axis xy equal tight;

        mp.F3.full.Nxi = size(mp.F3.full.mask.amp,2);
        mp.F3.full.Neta= size(mp.F3.full.mask.amp,1);   
        
        %--Number of points across the FPM in the compact model
        if(mp.F3.Rout==inf)
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
        
%         %--Coordinates for the FPMs in the full and compact models
%         if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.full.Nxi,2)==1  )
%             mp.F3.full.xisDL  = (-(mp.F3.full.Nxi -1)/2:(mp.F3.full.Nxi -1)/2)/mp.F3.full.res;
%             mp.F3.full.etasDL = (-(mp.F3.full.Neta-1)/2:(mp.F3.full.Neta-1)/2)/mp.F3.full.res;
%             
%             mp.F3.compact.xisDL  = (-(mp.F3.compact.Nxi -1)/2:(mp.F3.compact.Nxi -1)/2)/mp.F3.compact.res;
%             mp.F3.compact.etasDL = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2)/mp.F3.compact.res;
%         else
%             mp.F3.full.xisDL  = (-mp.F3.full.Nxi/2:(mp.F3.full.Nxi/2-1))/mp.F3.full.res;
%             mp.F3.full.etasDL = (-mp.F3.full.Neta/2:(mp.F3.full.Neta/2-1))/mp.F3.full.res;
%             
%             mp.F3.compact.xisDL  = (-mp.F3.compact.Nxi/2:(mp.F3.compact.Nxi/2-1))/mp.F3.compact.res;
%             mp.F3.compact.etasDL = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1))/mp.F3.compact.res;
%         end
        
        
        %--Make or read in focal plane mask (FPM) amplitude for the compact model
        FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
        mp.F3.compact.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
        
%         %--Downsample the FPM for the compact model --> DO NOT USE: Does not handle the true width of the mask correctly
%         [XIS,ETAS] = meshgrid(mp.F3.full.xisDL,mp.F3.full.etasDL);
%         [XIScompact,ETAScompact] = meshgrid(mp.F3.compact.xisDL,mp.F3.compact.etasDL);
%         if(mp.F3.Rout==inf); extrapval = 1; else extrapval = 0; end
%         mp.F3.compact.mask.amp = interp2(XIS,ETAS,mp.F3.full.mask.amp,XIScompact,ETAScompact,'cubic',extrapval);
%         % figure(205); imagesc(mp.F3.compact.mask.amp); axis xy equal tight;

end %--END OF FUNCTION