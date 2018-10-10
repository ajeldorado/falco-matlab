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

function mp = falco_config_gen_FPM_SPHLC(mp)

%%

        %--Generate iris part of the FPM at lower resolution.
        %  NOTE: The complex transmission value is multiplied within the SPHLC model.
        % Full model
        FPMgenInputs.pixresFPM = mp.F3.full.out.res; %--pixels per lambda_c/D
        FPMgenInputs.rhoInner = 0; % radius of inner FPM amplitude spot (in lambda_c/D)
        FPMgenInputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
        FPMgenInputs.FPMampFac = 1; % amplitude transmission of inner FPM spot
        FPMgenInputs.centering = mp.centering;
        mp.F3.full.iris = falco_gen_annular_FPM(FPMgenInputs);
        [mp.F3.full.out.Neta, mp.F3.full.out.Nxi] = size(mp.F3.full.iris);%--Number of points across the FPM in the full model
        % Compact model (some inputs the same as for full model)
        FPMgenInputs.pixresFPM = mp.F3.compact.out.res; %--pixels per lambda_c/D
        mp.F3.compact.iris = falco_gen_annular_FPM(FPMgenInputs);
        [mp.F3.compact.out.Neta, mp.F3.compact.out.Nxi] = size(mp.F3.compact.iris);%--Number of points across the FPM in the compact model

      
        %--Data for generating the complex-valued inner part of the FPM in the full and compact model functions:
        %--Number of points across the FPM in the full model
        mp.F3.full.in.Nxi  = mp.dm9.NdmPad;%ceil_even( 2*(mp.F3.Rin+buffer)*mp.F3.full.res);
        mp.F3.full.in.Neta = mp.F3.full.in.Nxi;
        %--Number of points across the FPM in the compact model
        mp.F3.compact.in.Nxi  = mp.dm9.compact.NdmPad; %ceil_even( 2*(mp.F3.Rin+buffer)*mp.F3.compact.res);
        mp.F3.compact.in.Neta = mp.F3.compact.in.Nxi;
        
        if(mp.flagSPHLConeFPM) %--Only inner region
            %--Generate iris part of the FPM at HIGHER resolution.  %  NOTE: The complex transmission value is multiplied within the SPHLC model.
            % Compact Model with 1-part propagation to/from F3 focal plane
            FPMgenInputs.pixresFPM = mp.F3.compact.in.res; %--pixels per lambda_c/D
            FPMgenInputs.rhoInner = 0; % radius of inner FPM amplitude spot (in lambda_c/D)
            FPMgenInputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
            FPMgenInputs.FPMampFac = 1; % amplitude transmission of inner FPM spot
            FPMgenInputs.centering = mp.centering;
            mp.F3.compact.irisHD = falco_gen_annular_FPM(FPMgenInputs);
            mp.F3.compact.irisHD = padOrCropEven(mp.F3.compact.irisHD,mp.F3.compact.in.Nxi);
            % Full Model with 1-part propagation to/from F3 focal plane
            FPMgenInputs.pixresFPM = mp.F3.full.in.res; %--pixels per lambda_c/D
            mp.F3.full.irisHD = falco_gen_annular_FPM(FPMgenInputs);
            mp.F3.full.irisHD = padOrCropEven(mp.F3.full.irisHD,mp.F3.full.in.Nxi);
        end
        
        
        %% Coordinates (dimensionless [DL]) for the FPMs in the full and compact models
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.full.in.Nxi,2)==1  )
            mp.F3.full.in.xisDL   = (-(mp.F3.full.in.Nxi -1)/2:(mp.F3.full.in.Nxi -1)/2)/mp.F3.full.in.res;
            mp.F3.full.out.xisDL  = (-(mp.F3.full.out.Nxi -1)/2:(mp.F3.full.out.Nxi -1)/2)/mp.F3.full.out.res;
            mp.F3.full.in.etasDL = mp.F3.full.in.xisDL;
            mp.F3.full.out.etasDL = mp.F3.full.out.xisDL;
            %mp.F3.full.etasDL = (-(mp.F3.full.Neta-1)/2:(mp.F3.full.Neta-1)/2)/mp.F3.full.res;
            
            mp.F3.compact.in.xisDL   = (-(mp.F3.compact.in.Nxi -1)/2:(mp.F3.compact.in.Nxi -1)/2)/mp.F3.compact.in.res;
            mp.F3.compact.out.xisDL  = (-(mp.F3.compact.out.Nxi -1)/2:(mp.F3.compact.out.Nxi -1)/2)/mp.F3.compact.out.res;
            mp.F3.compact.in.etasDL = mp.F3.compact.in.xisDL;
            mp.F3.compact.out.etasDL = mp.F3.compact.out.xisDL;
            %mp.F3.compact.etasDL = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2)/mp.F3.compact.res;
        else
            mp.F3.full.in.xisDL  = (-mp.F3.full.in.Nxi/2:(mp.F3.full.in.Nxi/2-1))/mp.F3.full.in.res;
            mp.F3.full.out.xisDL  = (-mp.F3.full.out.Nxi/2:(mp.F3.full.out.Nxi/2-1))/mp.F3.full.out.res;
            mp.F3.full.in.etasDL = mp.F3.full.in.xisDL;
            mp.F3.full.out.etasDL = mp.F3.full.out.xisDL;
            %mp.F3.full.etasDL = (-mp.F3.full.Neta/2:(mp.F3.full.Neta/2-1))/mp.F3.full.res;
            
            mp.F3.compact.in.xisDL  = (-mp.F3.compact.in.Nxi/2:(mp.F3.compact.in.Nxi/2-1))/mp.F3.compact.in.res;
            mp.F3.compact.out.xisDL  = (-mp.F3.compact.out.Nxi/2:(mp.F3.compact.out.Nxi/2-1))/mp.F3.compact.out.res;
            mp.F3.compact.in.etasDL = mp.F3.compact.in.xisDL;
            mp.F3.compact.out.etasDL = mp.F3.compact.out.xisDL;
            %mp.F3.compact.etasDL = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1))/mp.F3.compact.res;
        end
       


%% Coordinates (meters)

        mp.F3.full.in.dxi  = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.in.res;
        mp.F3.full.out.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.out.res;
        mp.F3.full.in.deta = mp.F3.full.in.dxi;
        mp.F3.full.out.deta = mp.F3.full.out.dxi;
        mp.F3.compact.in.dxi  = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.in.res;
        mp.F3.compact.out.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.out.res;
        mp.F3.compact.in.deta  = mp.F3.compact.in.dxi;
        mp.F3.compact.out.deta = mp.F3.compact.out.dxi;
        %--Compute coordinates in plane of FPM in the compact model (in meters)
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.compact.in.Nxi,2)==1  )
            %mp.xisFPMcompact = (-(mp.F3.compact.Nxi-1)/2:(mp.F3.compact.Nxi-1)/2)*mp.F3.compact.dxi;
            %mp.etasFPMcompact = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2).'*mp.F3.compact.deta;
            
            mp.F3.compact.in.xis   = (-(mp.F3.compact.in.Nxi -1)/2:(mp.F3.compact.in.Nxi -1)/2)*mp.F3.compact.in.dxi;
            mp.F3.compact.out.xis  = (-(mp.F3.compact.out.Nxi -1)/2:(mp.F3.compact.out.Nxi -1)/2)*mp.F3.compact.out.dxi;
            mp.F3.compact.in.etas  = (-(mp.F3.compact.in.Neta-1)/2:(mp.F3.compact.in.Neta-1)/2).'*mp.F3.compact.in.deta;
            mp.F3.compact.out.etas = (-(mp.F3.compact.out.Neta-1)/2:(mp.F3.compact.out.Neta-1)/2).'*mp.F3.compact.out.deta;
        else
            %mp.xisFPMcompact = (-mp.F3.compact.Nxi/2:(mp.F3.compact.Nxi/2-1))*mp.F3.compact.dxi;
            %mp.etasFPMcompact = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1)).'*mp.F3.compact.deta;
            
            mp.F3.compact.in.xis   = (-mp.F3.compact.in.Nxi/2: (mp.F3.compact.in.Nxi/2 -1))*mp.F3.compact.in.dxi;
            mp.F3.compact.out.xis  = (-mp.F3.compact.out.Nxi/2: (mp.F3.compact.out.Nxi/2 -1))*mp.F3.compact.out.dxi;
            mp.F3.compact.in.etas  = (-mp.F3.compact.in.Neta/2:(mp.F3.compact.in.Neta/2-1)).'*mp.F3.compact.in.deta;
            mp.F3.compact.out.etas = (-mp.F3.compact.out.Neta/2:(mp.F3.compact.out.Neta/2-1)).'*mp.F3.compact.out.deta;
        end
        
%         %--Focal Plane Coordinates for the SPHLC (in meters). Used for the DM9 Jacobian calculation.
%         if(strcmpi(mp.centering,'interpixel'))
%             mp.xisFPMinnerCompactCrop  = ( -(mp.NxiFPMinnerCompactCrop-1)/2:(mp.NxiFPMinnerCompactCrop-1)/2 )*mp.dxiFPMinnerCompact;
%             mp.etasFPMinnerCompactCrop  = ( -(mp.NetaFPMinnerCompactCrop-1)/2:(mp.NetaFPMinnerCompactCrop-1)/2 ).'*mp.detaFPMinnerCompact;        
%         else
%             mp.xisFPMinnerCompactCrop  = ( -mp.NxiFPMinnerCompactCrop/2:(mp.NxiFPMinnerCompactCrop/2-1) )*mp.dxiFPMinnerCompact;
%             mp.etasFPMinnerCompactCrop  = ( -mp.NetaFPMinnerCompactCrop/2:(mp.NetaFPMinnerCompactCrop/2-1) )*mp.detaFPMinnerCompact;
%         end
%         %--FPM Coordinates for SPHLC (in meters). Used for DM1 & DM2 Jacobians
%         mp.xisFPMinnerCompactMet = mp.xisFPMinnerCompact*(mp.fl*mp.lambda0/mp.P2.D);
%         mp.etasFPMinnerCompactMet = mp.etasFPMinnerCompact*(mp.fl*mp.lambda0/mp.P2.D);
%         mp.xisFPMouterCompactMet = mp.xisFPMouterCompact*(mp.fl*mp.lambda0/mp.P2.D);
%         mp.etasFPMouterCompactMet = mp.etasFPMouterCompact*(mp.fl*mp.lambda0/mp.P2.D);
%         
%         %--For plotting (in lambda_c/D)
%         mp.F3.compact.xisDL  = mp.xisFPMinnerCompactCrop/(mp.fl*mp.lambda0/mp.P2.D);
%         mp.F3.compact.etasDL = mp.etasFPMinnerCompactCrop/(mp.fl*mp.lambda0/mp.P2.D);

%%
end %--END OF FUNCTION
