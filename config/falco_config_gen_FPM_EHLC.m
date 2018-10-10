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

function mp = falco_config_gen_FPM_EHLC(mp)


        %% Generate Outer Iris part of the mask
        
        %--Generate outer opaque iris for FPM for the full model
        FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
        FPMgenInputs.rhoInner = 0; % radius of inner FPM amplitude spot (in lambda_c/D)
        FPMgenInputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
        FPMgenInputs.FPMampFac = 1; % amplitude transmission of inner FPM spot
        FPMgenInputs.centering = mp.centering;
        mp.F3.full.iris = falco_gen_annular_FPM(FPMgenInputs);
        %figure(204); imagesc(mp.F3.full.iris); axis xy equal tight; 
        
       %--Generate outer opaque iris for FPM for the compact model
        FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
        mp.F3.compact.iris = falco_gen_annular_FPM(FPMgenInputs);
        %figure(205); imagesc(mp.F3.compact.iris); axis xy equal tight; 

        %--Number of points across the total FPM (including the outer iris)
        mp.F3.full.Nxi = size(mp.F3.full.iris,2);
        mp.F3.full.Neta= size(mp.F3.full.iris,1);  
        mp.F3.compact.Nxi = size(mp.F3.compact.iris,2);
        mp.F3.compact.Neta= size(mp.F3.compact.iris,1);  
        
     
        %% PMGI bias layer (not necessarily uniform)
        
        DM8surfFull = falco_gen_EHLC_FPM_surf_from_cube(mp.dm8,'full');
        DM9surfFull = falco_gen_EHLC_FPM_surf_from_cube(mp.dm9,'full');

        if(mp.flagPlot); figure(8); imagesc(DM8surfFull); axis xy equal tight; colorbar; title('FPM Metal Thickness','Fontsize',20); set(gca,'Fontsize',20); drawnow; end
        if(mp.flagPlot); figure(9); imagesc(DM9surfFull); axis xy equal tight; colorbar; title('FPM Dielectric Thickness','Fontsize',20); set(gca,'Fontsize',20); drawnow; end
    
        DM8surfCompact = falco_gen_EHLC_FPM_surf_from_cube(mp.dm8,'compact');
        %DM9surfCompact = falco_gen_EHLC_FPM_surf_from_cube(mp.dm9,'compact');
        
        %--Dielectric (DM9) offsets
        switch mp.F3.dielProfile
            case{'uniform'}
                mp.F3.full.t_diel_nom_nm = mp.t_diel_bias_nm*ones(mp.F3.full.Nxi);
                mp.F3.compact.t_diel_nom_nm = mp.t_diel_bias_nm*ones(mp.F3.compact.Nxi);

            case{'flattop'}
                mp.F3.full.t_diel_nom_nm = mp.t_diel_bias_nm*ones(mp.F3.full.Nxi) - padOrCropEven(DM8surfFull*1e9,mp.F3.full.Nxi);
                mp.F3.compact.t_diel_nom_nm = mp.t_diel_bias_nm*ones(mp.F3.compact.Nxi) - padOrCropEven(DM8surfCompact*1e9,mp.F3.compact.Nxi);
                
            case{'pistondisk'}
                
                %--Make or read in DM8 disk for the full model
                FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
                FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
                FPMgenInputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
                FPMgenInputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
                FPMgenInputs.centering = mp.centering;
                diskFull =  1-falco_gen_annular_FPM(FPMgenInputs);
                %--Make or read in DM8 disk for the compact model
                FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
                diskCompact =  1-falco_gen_annular_FPM(FPMgenInputs);
                
                mp.F3.full.t_diel_nom_nm = mp.t_diel_bias_nm*ones(mp.F3.full.Nxi) + mp.F3.dielDiskPiston*padOrCropEven(diskFull,mp.F3.full.Nxi);
                mp.F3.compact.t_diel_nom_nm = mp.t_diel_bias_nm*ones(mp.F3.compact.Nxi) + mp.F3.dielDiskPiston*padOrCropEven(diskCompact,mp.F3.compact.Nxi);
                

                
            case{'ripple'}
                
                switch mp.centering
                    case 'pixel'
                        xsFull = (-mp.F3.full.Nxi/2:mp.F3.full.Nxi/2-1)/mp.F3.full.Nxi*2*mp.F3.Rout;
                        xsCompact = (-mp.F3.compact.Nxi/2:mp.F3.compact.Nxi/2-1)/mp.F3.compact.Nxi*2*mp.F3.Rout;
                    case 'interpixel'
                        xsFull = (-(mp.F3.full.Nxi-1)/2:(mp.F3.full.Nxi-1)/2)/mp.F3.full.Nxi*2*mp.F3.Rout;
                        xsCompact = (-(mp.F3.compact.Nxi-1)/2:(mp.F3.compact.Nxi-1)/2)/mp.F3.compact.Nxi*2*mp.F3.Rout;
                end
                
                [XSfull,YSfull] = meshgrid(xsFull);
                THETASfull = atan2(YSfull,XSfull);
                RSfull = sqrt(XSfull.^2 + YSfull.^2);
                %THETASfull(RSfull <= mp.F3.RminRipple) = 0;
                
                [XScompact,YScompact] = meshgrid(xsCompact);
                THETAScompact = atan2(YScompact,XScompact);
                RScompact = sqrt(XScompact.^2 + YScompact.^2);
                %THETAScompact(RScompact <= mp.F3.RminRipple) = 0;
                
                rippleFull = mp.F3.rippleAmp*cos( (THETASfull + deg2rad(mp.F3.ripplePhaseDeg))*mp.F3.rippleCharge);
                rippleFull(RSfull <= mp.F3.RminRipple) = 0;
                mp.F3.full.t_diel_nom_nm = rippleFull + mp.t_diel_bias_nm - padOrCropEven(DM8surfFull*1e9,mp.F3.full.Nxi);
                
                rippleCompact = mp.F3.rippleAmp*cos((THETAScompact+deg2rad(mp.F3.ripplePhaseDeg))*mp.F3.rippleCharge);
                rippleCompact(RScompact <= mp.F3.RminRipple) = 0;
                mp.F3.compact.t_diel_nom_nm = rippleCompact + mp.t_diel_bias_nm - padOrCropEven(DM8surfCompact*1e9,mp.F3.compact.Nxi);
  
        end
        if(mp.flagPlot); figure(19); imagesc(mp.F3.full.t_diel_nom_nm); axis xy equal tight; colorbar; title('Starting FPM Dielectric','Fontsize',20); set(gca,'Fontsize',20); drawnow; end
        if(mp.flagPlot); figure(20); imagesc(mp.F3.compact.t_diel_nom_nm); axis xy equal tight; colorbar; title('Starting FPM Dielectric','Fontsize',20); set(gca,'Fontsize',20); drawnow; end
        


end %--END OF FUNCTION
