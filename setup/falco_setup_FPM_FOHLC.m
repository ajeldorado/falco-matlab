% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%
% REVISION HISTORY:
% --------------
% Modified on 2018-03-22 by A.J. Riggs from HLC to FOHLC.
% Created on 2018-10-01 by A.J. Riggs by extracting material from falco_init_ws.m.
% ---------------

function mp = falco_setup_FPM_FOHLC(mp)

%% DM8 and DM9 (Optimizable FPM) Setup

%--Centering of DM surfaces on array
mp.dm8.centering = mp.centering;
mp.dm9.centering = mp.centering;

mp.dm9.compact = mp.dm9;

%%
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
        
%% Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)
r_cent_lam0D = mp.dm9.r_cent_act*mp.dm9.dm_spacing/(mp.dm9.dxi)/mp.F3.full.res;

%%
mp.F3.RinA_inds = [];
mp.F3.RinAB_inds = [];
for ii=1:mp.dm9.NactTotal
    %--Zero out FPM actuators beyond the allowed radius (mp.F3.Rin)
    if(r_cent_lam0D(ii) > mp.F3.Rin-mp.dm9.FPMbuffer)
       mp.dm9.inf_datacube(:,:,ii) = 0*mp.dm9.inf_datacube(:,:,ii);
       mp.dm9.compact.inf_datacube(:,:,ii) = zeros(size(mp.dm9.compact.inf_datacube(:,:,ii)));
    end

    %--Get the indices for the actuators within radius mp.F3.RinA
    if(r_cent_lam0D(ii) <= mp.F3.RinA-mp.dm9.FPMbuffer)
        mp.F3.RinA_inds = [mp.F3.RinA_inds; ii];
    else %--Get the indices for the actuators between radii mp.F3.RinA and mp.F3.Rin
        mp.F3.RinAB_inds = [mp.F3.RinAB_inds; ii];
    end
end
fprintf('%d actuators in DM9.\n',mp.dm9.NactTotal);

%--Uniform gain for all actuators
mp.dm9.VtoHavg = 1e-9; %--gain of DM9 (meters of phase/Volt)
mp.dm9.VtoH = mp.dm9.VtoHavg*ones(mp.dm9.NactTotal,1); % Gains: volts to meters in surface height;

if(isfield(mp.dm9,'V')==false); mp.dm9.V = zeros(mp.dm9.NactTotal,1); mp.dm9.V(mp.F3.RinA_inds) = mp.dm9.V0coef; else; mp.dm9.V = mp.DM9V0; end %--Initial DM9 voltages

%--DM8 Option 1: Set as same basis set (DM actuators) as DM9. BE CAREFUL ABOUT OVERWRITING VALUES INADVERTENTLY!!!!!!!!!!!
dm8Vmin = mp.dm8.Vmin;
dm8Vmax = mp.dm8.Vmax;
dm8v0coef = mp.dm8.V0coef; %--Starting value. Save as temporary value to avoid overwriting
dm8VtoHavg = mp.dm8.VtoHavg; %--Starting value. Save as temporary value to avoid overwriting
dm8actSens = mp.dm8.act_sens;
dm8weight = mp.dm8.weight;
dm8maxAbsdV = mp.dm8.maxAbsdV;%--Starting value. Save as temporary value to avoid overwriting
mp.dm8 = mp.dm9; %--DANGER: OVERWRITING ANY PREVIOUS mp.dm8 DATA!!! -----------------------------------------------!!!!!!!!!!!!!!!!!!!!
mp.dm8 = rmfield(mp.dm8,'V');
mp.dm8.V0coef = dm8v0coef;
mp.dm8.VtoHavg = dm8VtoHavg;
mp.dm8.Vmin = dm8Vmin;
mp.dm8.Vmax = dm8Vmax;
mp.dm8.act_sens = dm8actSens;
mp.dm8.weight = dm8weight;
mp.dm8.maxAbsdV = dm8maxAbsdV;%--Starting value. Save as temporary value to avoid overwriting
fprintf('%d actuators in DM8.\n',mp.dm8.NactTotal);
mp.dm8.VtoH = mp.dm8.VtoHavg*ones(mp.dm8.NactTotal,1); % Gains: volts to meters in surface height;
if(isfield(mp,'DM8V0')) %--Initial DM voltages
    mp.dm8.V = mp.DM8V0;
else
    mp.dm8.V = zeros(mp.dm8.NactTotal,1);
    mp.dm8.V(mp.F3.RinA_inds) = mp.dm8.V0coef; 
end

% %--Zero out parts of DM9 actuators that go outside the nickel disk. Also apply the grayscale edge.
for ii=1:mp.dm9.NactTotal
    if(r_cent_lam0D(ii) > mp.F3.RinA-mp.dm9.FPMbuffer)
        %--Zero out FPM actuators beyond the inner spot (mp.F3.Rin)
        mp.dm9.inf_datacube(:,:,ii) = 0*mp.dm9.inf_datacube(:,:,ii);
        mp.dm9.compact.inf_datacube(:,:,ii) = zeros(size(mp.dm9.compact.inf_datacube(:,:,ii)));
    end
end
fprintf('%d actuators in DM9.\n',mp.dm9.NactTotal);

end %--END OF FUNCTION