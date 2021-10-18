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

function mp = falco_setup_FPM_HLC(mp)

mp.dm9.Nact = ceil_even(2*mp.F3.Rin*mp.dm9.actres); % number of actuators across DM9 (if not in a hex grid)

switch mp.dm9.inf0name
    case '3x3'
        mp.dm9.inf0 = 1/4*[1, 2, 1; 2, 4, 2; 1, 2, 1]; % influence function
        mp.dm9.dx_inf0_act = 1/2; % number of inter-actuator widths per pixel 
        %--FPM resolution (pixels per lambda0/D) in the compact and full models.
        mp.F3.compact.res = mp.dm9.actres/mp.dm9.dx_inf0_act;
        mp.F3.full.res = mp.dm9.actres/mp.dm9.dx_inf0_act;

    case 'Lanczos3'
        N = 91;
        xs = (-(N-1)/2:(N-1)/2)/N*10*0.906;

        a = 3;
        Lx0 = a*sin(pi*xs).*sin(pi*xs/a)./(pi*xs).^2;
        Lx0(xs==0) = 1;
        Lx = Lx0;
        Lx(xs>=a) = 0;
        Lx(xs<=-a) = 0;
        Lxy = Lx.'*Lx; %--The 2-D Lanczos kernel
        Nhalf = ceil(N/2);
        Nrad = 30;
        inds = Nhalf-Nrad:Nhalf+Nrad;
        LxyCrop = Lxy(inds,inds);

        mp.dm9.inf0 = LxyCrop; % influence function
        mp.dm9.dx_inf0_act = 1/10; % number of inter-actuator widths per pixel 

    case 'Xinetics'
        mp.dm9.inf0 = 1*fitsread('influence_dm5v2.fits');
        mp.dm9.dx_inf0_act = 1/10; % number of inter-actuator widths per pixel 
end
    
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

mp.dm9.ABfac = 1; %--Gain factor between inner and outer FPM regions
mp.dm9.VtoHavg = 1e-9; %--gain of DM9 (meters/Volt)
mp.dm9.VtoH = mp.dm9.VtoHavg*ones(mp.dm9.NactTotal,1); % Gains: volts to meters in surface height;
mp.dm9.VtoH(mp.F3.RinAB_inds) = mp.dm9.ABfac*mp.dm9.VtoH(mp.F3.RinAB_inds);

if(isfield(mp.dm9,'V')==false) %--Initial DM9 voltages
    mp.dm9.V = zeros(mp.dm9.NactTotal,1);
    mp.dm9.V(mp.F3.RinA_inds) = mp.dm9.V0coef;
else
    mp.dm9.V = mp.DM9V0;
end

mp.dm9.Vmin = min(mp.t_diel_nm_vec); % minimum thickness of FPM dielectric layer (nm)
mp.dm9.Vmax = max(mp.t_diel_nm_vec); % maximum thickness (from one actuator, not of the facesheet) of FPM dielectric layer (nm)

%%%%%---OPTIONS FOR DEFINING DM8 (FPM Metal)
mp.dm8.VtoHavg = 1e-9; %--gain of DM8 (meters/Volt)
mp.dm8.Vmin = min(mp.t_metal_nm_vec); % minimum thickness of FPM metal layer (nm)
mp.dm8.Vmax = max(mp.t_metal_nm_vec); % maximum thickness (from one actuator, not of the facesheet) of FPM metal layer (nm)

%--DM8 Option 2: Set basis as a single nickel disk.
mp.dm8.NactTotal = 1;
mp.dm8.act_ele = 1;
fprintf('%d actuators in DM8.\n',mp.dm8.NactTotal);
mp.dm8.VtoH = mp.dm8.VtoHavg*ones(mp.dm8.NactTotal,1); % Gains: volts to meters in surface height;
mp.dm8.xy_box_lowerLeft = [1;1];
mp.dm8.compact = mp.dm8;
if(isfield(mp.dm8,'V')==false) %--Initial DM voltages
    mp.dm8.V = mp.dm8.V0coef*ones(mp.dm8.NactTotal,1);
else
    mp.dm8.V = mp.DM8V0;
end
% Don't define extra actuators and time:
if(mp.F3.Rin~=mp.F3.RinA); error('falco_init_ws.m: Change mp.F3.Rin and mp.F3.RinA to be equal to avoid wasting time.'); end

% Copy over some common values from DM9:
mp.dm8.dxi = mp.dm9.dxi; %--Width of a pixel at the FPM in full model (meters)
mp.dm8.NdmPad = mp.dm9.NdmPad;
mp.dm8.Nbox = mp.dm8.NdmPad;
mp.dm8.compact.dxi = mp.dm9.compact.dxi; %--Width of a pixel at the FPM in compact model (meters)
mp.dm8.compact.NdmPad = mp.dm9.compact.NdmPad;
mp.dm8.compact.Nbox = mp.dm8.compact.NdmPad;

%--Make or read in DM8 disk for the full model
FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
FPMgenInputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
FPMgenInputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
FPMgenInputs.centering = mp.centering;
mp.dm8.inf_datacube =  round(pad_crop(1-falco_gen_annular_FPM(FPMgenInputs), mp.dm8.NdmPad));
%--Make or read in DM8 disk for the compact model
FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
mp.dm8.compact.inf_datacube =  round(pad_crop(1-falco_gen_annular_FPM(FPMgenInputs), mp.dm8.compact.NdmPad));

%--Zero out parts of DM9 actuators that go outside the nickel disk. Also apply the grayscale edge.
DM8windowFull = mp.dm8.inf_datacube;
DM8windowCompact = mp.dm8.compact.inf_datacube;
for iact=1:mp.dm9.NactTotal 
    if( any(any(mp.dm9.inf_datacube(:,:,iact))) && any(mp.dm9.VtoH(iact)) )
        y_box_ind = mp.dm9.xy_box_lowerLeft(1,iact):mp.dm9.xy_box_lowerLeft(1,iact)+mp.dm9.Nbox-1; % x-indices in pupil arrays for the box
        x_box_ind = mp.dm9.xy_box_lowerLeft(2,iact):mp.dm9.xy_box_lowerLeft(2,iact)+mp.dm9.Nbox-1; % y-indices in pupil arrays for the box
        mp.dm9.inf_datacube(:,:,iact) = DM8windowFull(x_box_ind,y_box_ind).*mp.dm9.inf_datacube(:,:,iact);

        y_box_ind = mp.dm9.compact.xy_box_lowerLeft(1,iact):mp.dm9.compact.xy_box_lowerLeft(1,iact)+mp.dm9.compact.Nbox-1; % x-indices in pupil arrays for the box
        x_box_ind = mp.dm9.compact.xy_box_lowerLeft(2,iact):mp.dm9.compact.xy_box_lowerLeft(2,iact)+mp.dm9.compact.Nbox-1; % y-indices in pupil arrays for the box
        mp.dm9.compact.inf_datacube(:,:,iact) = DM8windowCompact(x_box_ind,y_box_ind).*mp.dm9.compact.inf_datacube(:,:,iact);
    end
end

end %--END OF FUNCTION