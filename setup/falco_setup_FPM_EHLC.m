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

function mp = falco_setup_FPM_EHLC(mp)


%--Centering of DM surfaces on array
mp.dm8.centering = mp.centering;
mp.dm9.centering = mp.centering;

% if(isfield(mp.dm9,'inf_datacube')==0 && any(mp.dm_ind==9) )
    %mp.dm9.dx_inf0_act = 1/10; % 'influence_dm5v2.fits' units of actuator spacings, sampling of the influence function. Always = 1/10 for 
mp.dm9.compact = mp.dm9;


%%

        %--DM9
        if(isfield(mp,'flagDM9inf3x3'))
            mp.dm9.xcent_dm = mp.dm9.Nact/2 - 1/2;
            mp.dm9.ycent_dm = mp.dm9.Nact/2 - 1/2;
            if(strcmpi(mp.centering,'interpixel'))
               error('falco_init_ws: The 3x3 influence function for DM9 requires a pixel-centered coordinate system (for now).');
            end
        else
            mp.dm9.xcent_dm = mp.dm9.Nact/2 - 1/2;
            mp.dm9.ycent_dm = mp.dm9.Nact/2 - 1/2;
        end
        mp.dm9.dm_spacing = 1/mp.dm9.actres*(mp.fl*mp.lambda0/mp.P2.D); % meters, pitch of DM actuators
        mp.dm9.compact = mp.dm9;

        mp.dm9.compact.dm_spacing = mp.dm9.dm_spacing; % meters, pitch of DM actuators

        mp.dm9.dx_inf0 = (mp.dm9.dx_inf0_act)*mp.dm9.dm_spacing; % meters, sampling of the influence function 
        mp.dm9.compact.dx_inf0 = (mp.dm9.compact.dx_inf0_act)*mp.dm9.compact.dm_spacing; % meters, sampling of the influence function 
        
        mp.dm9.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.res; % width of a pixel at the FPM in the full model (meters)
        mp.dm9.compact.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.res; % width of a pixel at the FPM in the compact model (meters)
        
        if(length(mp.dm9.inf0)==3)
            mp.dm9 = falco_fpm_inf_cube_3x3(mp.dm9);
            mp.dm9.compact = falco_fpm_inf_cube_3x3(mp.dm9.compact);
        else
            mp.dm9 = falco_fpm_inf_cube(mp.dm9);
            mp.dm9.compact = falco_fpm_inf_cube(mp.dm9.compact);
        end
        
        
        %%--DM8 (FPM Metal)
        
        mp.dm8.VtoHavg = 1e-9; %--gain of DM8 (meters/Volt)
        
        if( strcmpi(mp.F3.metalProfile,'tophat') && mp.F3.flagFixedDisk==true )
            
            %--DM8 Option 2: Set basis as a single nickel disk.
            mp.dm8.NactTotal = 1;
            mp.dm8.act_ele = 1;
            fprintf('%d actuators in DM8.\n',mp.dm8.NactTotal);
            mp.dm8.VtoH = mp.dm8.VtoHavg*ones(mp.dm8.NactTotal,1); % Gains: volts to meters in surface height;
            mp.dm8.xy_box_lowerLeft = [1;1];
            mp.dm8.compact = mp.dm8;
            if(isfield(mp.dm8,'V')==false); mp.dm8.V = mp.F3.metalPeak*ones(mp.dm8.NactTotal,1); else; mp.dm8.V = DM8V0; end %--Initial DM voltages
% %             if(isfield(mp.dm8,'V')==false); mp.dm8.V = mp.dm8.V0coef*ones(mp.dm8.NactTotal,1); else; mp.dm8.V = DM8V0; end %--Initial DM voltages
            % Don't define extra actuators and time:
            %if(mp.F3.Rin~=mp.F3.RinA); error('falco_init_ws.m: Change mp.F3.Rin and mp.F3.RinA to be equal to avoid wasting time.'); end

            % Copy over some common values from DM9:
            mp.dm8.dxi = mp.dm9.dxi; %--Width of a pixel at the FPM in full model (meters)
            mp.dm8.compact.dxi = mp.dm9.compact.dxi; %--Width of a pixel at the FPM in compact model (meters)

            %--Make or read in DM8 disk for the full model
            FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
            FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
            FPMgenInputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
            FPMgenInputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
            FPMgenInputs.centering = mp.centering;
            mp.dm8.inf_datacube =  1-falco_gen_annular_FPM(FPMgenInputs);
            mp.dm8.NdmPad = length(mp.dm8.inf_datacube);
            mp.dm8.Nbox = mp.dm8.NdmPad;
            %--Make or read in DM8 disk for the compact model
            FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
            mp.dm8.compact.inf_datacube =  1-falco_gen_annular_FPM(FPMgenInputs);
            mp.dm8.compact.NdmPad = length(mp.dm8.compact.inf_datacube);
            mp.dm8.compact.Nbox = mp.dm8.compact.NdmPad;
            % figure(204); imagesc(mp.dm8.inf_datacube); axis xy equal tight; colorbar; drawnow;
            % figure(205); imagesc(mp.dm8.compact.inf_datacube); axis xy equal tight; colorbar; drawnow;

            
            
        else
            
            %--DM8 (same influence function and sampling, but can have different radius/number of actuators
            mp.dm8.xcent_dm = mp.dm8.Nact/2 - 1/2; %mp.dm9.xcent_dm;
            mp.dm8.ycent_dm = mp.dm8.Nact/2 - 1/2;%mp.dm9.ycent_dm;

            mp.dm8.dm_spacing = mp.dm9.dm_spacing;
            mp.dm8.inf0 = mp.dm9.inf0;

            mp.dm8.compact = mp.dm8;

            mp.dm8.dx_inf0_act = mp.dm9.dx_inf0_act;
            mp.dm8.compact.inf0 = mp.dm9.compact.inf0;
            mp.dm8.compact.dx_inf0_act = mp.dm9.compact.dx_inf0_act;

            mp.dm8.dx_inf0 = mp.dm9.dx_inf0;
            mp.dm8.compact.dx_inf0 = mp.dm9.compact.dx_inf0;

            mp.dm8.dxi = mp.dm9.dxi ;
            mp.dm8.compact.dxi = mp.dm9.compact.dxi;

            if(length(mp.dm8.inf0)==3)
                mp.dm8 = falco_fpm_inf_cube_3x3(mp.dm8);
                mp.dm8.compact = falco_fpm_inf_cube_3x3(mp.dm8.compact);
            else
                mp.dm8 = falco_fpm_inf_cube(mp.dm8);
                mp.dm8.compact = falco_fpm_inf_cube(mp.dm8.compact);
            end
        
            %%%%%---OPTIONS FOR DEFINING DM8 (FPM Metal)
            mp.dm8.VtoH = mp.dm8.VtoHavg*ones(mp.dm8.NactTotal,1);%1*1e-9*ones(mp.dm9.Nact); % Gains: volts to meters in surface height;

            mp.dm8.Vmin = min(mp.t_metal_nm_vec); % minimum thickness of FPM metal layer (nm)
            mp.dm8.Vmax = max(mp.t_metal_nm_vec); % maximum thickness (from one actuator, not of the facesheet) of FPM metal layer (nm)

            %--Initial DM8 voltages
            % Approximate conversion from amplitude to Ni thickness at 575nm:  t_Ni = log10(amp)*(183e-9 - 50e-9)/(-3 +0.861) 

            if(isfield(mp.dm8,'V')==false)
                mp.dm8.V = zeros(mp.dm8.NactTotal,1); 
                rs = mp.dm8.r_cent_act(:)*mp.dm8.dm_spacing/(mp.dm8.dxi)/mp.F3.full.res;

                switch lower(mp.F3.metalProfile)
                    case 'sinc4' 
                        mp.F3.b = 6.5;  % amplitude = (1-sinc(pi*rs/2/b).^2).^2;                    
                        ampVec = (1-sinc(pi*rs(:)/2/mp.F3.b).^2).^2;
                        mp.dm8.V = 1e9*log10(ampVec)*(183e-9 - 50e-9)/(-3 +0.861);

                    case{'tophat'}
                        mp.dm8.V(rs <= mp.F3.Rin) = mp.F3.metalPeak;
                end

            else
                mp.dm8.V = DM8V0; 
            end 

            %--Zero out influence functions for actuators outside a given radius
            r_cent_lam0D = mp.dm8.r_cent_act*mp.dm8.dm_spacing/(mp.dm8.dxi)/mp.F3.full.res;
            for ii=1:mp.dm8.NactTotal
                %--Zero out FPM actuators beyond the allowed radius (mp.F3.RinMaxDiel)
                if(r_cent_lam0D(ii) > mp.F3.RinMaxMetal) %-mp.dm8.FPMbuffer)
                   mp.dm8.inf_datacube(:,:,ii) = 0*mp.dm8.inf_datacube(:,:,ii);
                   mp.dm8.compact.inf_datacube(:,:,ii) = zeros(size(mp.dm8.compact.inf_datacube(:,:,ii)));
                end
            end

        end
        
        %--Coordinates for the DM8 sub-arrays (actuated regions) in the full and compact models
        if( strcmpi(mp.centering,'interpixel') )
            mp.dm8.xisDL = (-(mp.dm8.NdmPad-1)/2:(mp.dm8.NdmPad-1)/2)/mp.F3.full.res;
            mp.dm8.compact.xisDL = (-(mp.dm8.compact.NdmPad-1)/2:(mp.dm8.compact.NdmPad-1)/2)/mp.F3.compact.res;
        else
            mp.dm8.xisDL  = (-mp.dm8.NdmPad/2:(mp.dm8.NdmPad/2-1))/mp.F3.full.res;
            mp.dm8.compact.xisDL  = (-mp.dm8.compact.NdmPad/2:(mp.dm8.compact.NdmPad/2-1))/mp.F3.compact.res;
        end
        
        
        %--Coordinates for the DM9 sub-arrays (actuated regions) in the full and compact models
        if( strcmpi(mp.centering,'interpixel') )
            mp.dm9.xisDL = (-(mp.dm9.NdmPad-1)/2:(mp.dm9.NdmPad-1)/2)/mp.F3.full.res;
            mp.dm9.compact.xisDL = (-(mp.dm9.compact.NdmPad-1)/2:(mp.dm9.compact.NdmPad-1)/2)/mp.F3.compact.res;
        else
            mp.dm9.xisDL  = (-mp.dm9.NdmPad/2:(mp.dm9.NdmPad/2-1))/mp.F3.full.res;
            mp.dm9.compact.xisDL  = (-mp.dm9.compact.NdmPad/2:(mp.dm9.compact.NdmPad/2-1))/mp.F3.compact.res;
        end
        
        fprintf('%d actuators in DM8.\n',mp.dm8.NactTotal);
        fprintf('%d actuators in DM9.\n',mp.dm9.NactTotal);

        
        mp.dm9.VtoHavg = 1e-9; %--gain of DM9 (meters/Volt)
        mp.dm9.VtoH = mp.dm9.VtoHavg*ones(mp.dm9.NactTotal,1);%1*1e-9*ones(mp.dm9.Nact); % Gains: volts to meters in surface height;
%         mp.dm9.VtoH(mp.F3.RinAB_inds) = mp.dm9.ABfac*mp.dm9.VtoH(mp.F3.RinAB_inds);
%         if(isfield(mp.dm9,'V')==false); mp.dm9.V = zeros(mp.dm9.NactTotal,1); mp.dm9.V(mp.F3.RinA_inds) = mp.dm9.V0coef; else; mp.dm9.V = DM9V0; end %--Initial DM9 voltages
        
        mp.dm9.Vmin = min(mp.t_diel_nm_vec)-mp.t_diel_bias_nm;  % minimum thickness of FPM dielectric layer (nm)
        mp.dm9.Vmax = max(mp.t_diel_nm_vec);  % maximum thickness (from one actuator, not of the facesheet) of FPM dielectric layer (nm)
           
        %--DM9 voltages from offset/bias
        if(isfield(mp.dm9,'V')==false); mp.dm9.V = zeros(mp.dm9.NactTotal,1); else; mp.dm9.V = DM9V0; end %--Initial DM9 voltages

        %--Zero out influence functions for actuators outside a given radius
        r_cent_lam0D = mp.dm9.r_cent_act*mp.dm9.dm_spacing/(mp.dm9.dxi)/mp.F3.full.res;
        for ii=1:mp.dm9.NactTotal
            %--Zero out FPM actuators beyond the allowed radius (mp.F3.RinMaxDiel)
            if(r_cent_lam0D(ii) > mp.F3.RinMaxDiel) %-mp.dm9.FPMbuffer)
               mp.dm9.inf_datacube(:,:,ii) = 0*mp.dm9.inf_datacube(:,:,ii);
               mp.dm9.compact.inf_datacube(:,:,ii) = zeros(size(mp.dm9.compact.inf_datacube(:,:,ii)));
            end
        end







%%
end %--END OF FUNCTION
