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

function mp = falco_setup_FPM_HLC_3foldZern(mp)


%% DM8 and DM9 (Optimizable FPM) Setup


%--Centering of DM surfaces on array
mp.dm8.centering = mp.centering;
mp.dm9.centering = mp.centering;

mp.dm9.compact = mp.dm9;

% mp.dm9.FPMbuffer = 0; %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)

%% DM9 as Zernikes

%--Pixel size [meters]
mp.dm9.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.res; % width of a pixel at the FPM in the full model (meters)
mp.dm9.compact.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.res; % width of a pixel at the FPM in the compact model (meters)
 

% if(isfield(mp,'flagDM9inf3x3'))
%     mp.dm9.xcent_dm = mp.dm9.Nact/2 - 1/2;
%     mp.dm9.ycent_dm = mp.dm9.Nact/2 - 1/2;
%     if(strcmpi(mp.centering,'interpixel'))
%        error('falco_init_ws: The 3x3 influence function for DM9 requires a pixel-centered coordinate system (for now).');
%     end
% else
%     mp.dm9.xcent_dm = mp.dm9.Nact/2 - 1/2;
%     mp.dm9.ycent_dm = mp.dm9.Nact/2 - 1/2;
% end
% mp.dm9.dm_spacing = 1/mp.dm9.actres*(mp.fl*mp.lambda0/mp.P2.D); % meters, pitch of DM actuators
% mp.dm9.compact = mp.dm9;

% mp.dm9.compact.dm_spacing = mp.dm9.dm_spacing; % meters, pitch of DM actuators

% mp.dm9.dx_inf0 = (mp.dm9.dx_inf0_act)*mp.dm9.dm_spacing;%1e-4; % meters, sampling of the influence function 

% mp.dm9.compact.dx_inf0 = (mp.dm9.compact.dx_inf0_act)*mp.dm9.compact.dm_spacing;%1e-4; % meters, sampling of the influence function 


% if(length(mp.dm9.inf0)==3)
%     mp.dm9 = falco_fpm_inf_cube_3x3(mp.dm9);
%     mp.dm9.compact = falco_fpm_inf_cube_3x3(mp.dm9.compact);
% else
%     mp.dm9 = falco_fpm_inf_cube(mp.dm9);
%     mp.dm9.compact = falco_fpm_inf_cube(mp.dm9.compact);
% end

%--Generate datacube of influence functions, which are Zernikes
NbeamCompact = mp.F3.compact.res*mp.F3.Rin*2;
mp.dm9.compact.inf_datacube = falco_3fold_symmetry_Zernikes(NbeamCompact,mp.dm9.maxRadialOrder,mp.centering,'SymmAxis','y');
mp.dm9.compact.NdmPad = size(mp.dm9.compact.inf_datacube,1);
mp.dm9.compact.Nbox = mp.dm9.compact.NdmPad; %--the Zernikes take up the full array.

NbeamFull = mp.F3.full.res*mp.F3.Rin*2;
mp.dm9.inf_datacube = falco_3fold_symmetry_Zernikes(NbeamFull,mp.dm9.maxRadialOrder,mp.centering,'SymmAxis','y');
mp.dm9.NdmPad = size(mp.dm9.inf_datacube,1);
mp.dm9.Nbox = mp.dm9.NdmPad; %--the Zernikes take up the full array.

mp.dm9.NactTotal = size(mp.dm9.inf_datacube,3);
mp.dm9.VtoH = mp.dm9.VtoHavg*ones(mp.dm9.NactTotal,1);

%--Lower-left pixel coordinates are all (1,1) since the Zernikes take up the full array.
mp.dm9.xy_box_lowerLeft = ones(2,mp.dm9.NactTotal);
mp.dm9.compact.xy_box_lowerLeft = ones(2,mp.dm9.NactTotal);


%--Coordinates for the full FPM array
if(strcmpi(mp.centering,'pixel')  ) 
    mp.dm9.compact.x_pupPad = (-mp.dm9.compact.NdmPad/2:(mp.dm9.compact.NdmPad/2 - 1))*mp.dm9.compact.dxi; % meters, coords for the full DM arrays. Origin is centered on a pixel
else
    mp.dm9.compact.x_pupPad = (-(mp.dm9.compact.NdmPad-1)/2:(mp.dm9.compact.NdmPad-1)/2)*mp.dm9.compact.dxi; % meters, coords for the full DM arrays. Origin is centered between pixels for an even-sized array
end
mp.dm9.compact.y_pupPad = mp.dm9.compact.x_pupPad;

%%
%--Initial DM9 voltages
if(isfield(mp.dm9,'V')==false)
    inf1 = mp.dm9.inf_datacube(:,:,1);
    meanVal = mean(inf1(inf1~=0));
    mp.dm9.V = zeros(mp.dm9.NactTotal,1); 
    mp.dm9.V(1) = mp.dm9.V0coef/meanVal; 
else
    mp.dm9.V = mp.DM9V0; 
end 

mp.dm9.Vmin = min(mp.t_diel_nm_vec);  % minimum thickness of FPM dielectric layer (nm)
mp.dm9.Vmax = max(mp.t_diel_nm_vec);  % maximum thickness (from one actuator, not of the facesheet) of FPM dielectric layer (nm)

%% Apply window to DM9 influence functions

hg_expon = 44; %--Found empirically
apRad = mp.F3.Rin/(mp.F3.Rin+0.1); %--Found empirically
OD = 1;

%--Full Model
Narray = mp.dm9.NdmPad;
Nbeam = NbeamFull;
%--Coordinates normalized to the beam radius (not diameter)
switch mp.centering
    case 'interpixel'
        xs = (-(Narray-1)/2:(Narray-1)/2)/Nbeam*2;
    case 'pixel'
        xs = (-Narray/2:(Narray/2-1))/Nbeam*2;
end
[XS,YS] = meshgrid(xs);
RS = sqrt(XS.^2+YS.^2); 
THETAS = atan2(YS,XS);
mask = RS<=1;
windowFull = mask.*exp(-(RS/(apRad*OD)).^hg_expon);

% figure(21); imagesc(mask); axis xy equal tight; colorbar;
% figure(22); imagesc(windowFull); axis xy equal tight; colorbar;
% figure(23); plot(xs*mp.F3.Rin,windowFull(Narray/2,:)); 
% ylim([-0.05,1.05])
% xlim(mp.F3.Rin*[0.87,1.1])

%--Compact Model
Narray = mp.dm9.compact.NdmPad;
Nbeam = NbeamCompact;
%--Coordinates normalized to the beam radius (not diameter)
switch mp.centering
    case 'interpixel'
        xs = (-(Narray-1)/2:(Narray-1)/2)/Nbeam*2;
    case 'pixel'
        xs = (-Narray/2:(Narray/2-1))/Nbeam*2;
end
[XS,YS] = meshgrid(xs);
RS = sqrt(XS.^2+YS.^2); 
THETAS = atan2(YS,XS);
mask = RS<=1;
windowCompact = mask.*exp(-(RS/(apRad*OD)).^hg_expon);


for di = 1:mp.dm9.NactTotal
%     figure(41); imagesc(mp.dm9.inf_datacube(:,:,di)); axis xy equal tight; colorbar; drawnow;
    mp.dm9.inf_datacube(:,:,di) = windowFull.*mp.dm9.inf_datacube(:,:,di);
    mp.dm9.compact.inf_datacube(:,:,di) = windowCompact.*mp.dm9.compact.inf_datacube(:,:,di);
%     figure(42); imagesc(mp.dm9.inf_datacube(:,:,di)); axis xy equal tight; colorbar; drawnow;
%     pause(0.1);
end

%% DM8 as a disk (piston as only influence function)

%%%%%---OPTIONS FOR DEFINING DM8 (FPM Metal)
mp.dm8.VtoHavg = 1e-9; %--gain of DM8 (meters/Volt) %--MOVE TO overwritable defaults when controlling DM8.
mp.dm8.Vmin = min(mp.t_metal_nm_vec); % minimum thickness of FPM metal layer (nm)
mp.dm8.Vmax = max(mp.t_metal_nm_vec); % maximum thickness (from one actuator, not of the facesheet) of FPM metal layer (nm)

%         %--DM8 Option 1: Set as same basis set (DM actuators) as DM9. 
%         dm8v0coef = mp.dm8.V0coef; %--Starting value. Save as temporary value to avoid overwriting
%         mp.dm8 = mp.dm9; %--DANGER: OVERWRITING ANY PREVIOUS mp.dm8 DATA!!!
%         mp.dm8.V0coef = dm8v0coef;
%         fprintf('%d actuators in DM8.\n',mp.dm8.NactTotal);
%         mp.dm8.ABfac = 1;%1/5; %--Gain factor between inner and outer FPM regions
%         mp.dm8.VtoH = mp.dm8.VtoHavg*ones(mp.dm8.NactTotal,1);%1*1e-9*ones(mp.dm9.Nact); % Gains: volts to meters in surface height;
%         mp.dm8.VtoH(mp.F3.RinAB_inds) = mp.dm8.ABfac*mp.dm8.VtoH(mp.F3.RinAB_inds);
%         if(isfield(mp.dm8,'V')==false); mp.dm8.V = zeros(mp.dm8.NactTotal,1); mp.dm8.V(mp.F3.RinA_inds) = mp.dm8.V0coef; else; mp.dm8.V = DM8V0; end %--Initial DM voltages


%--DM8 Option 2: Set basis as a single nickel disk.
mp.dm8.NactTotal = 1;
mp.dm8.act_ele = 1;
fprintf('%d actuators in DM8.\n',mp.dm8.NactTotal);
mp.dm8.VtoH = mp.dm8.VtoHavg*ones(mp.dm8.NactTotal,1); % Gains: volts to meters in surface height;
mp.dm8.xy_box_lowerLeft = [1;1];
mp.dm8.compact = mp.dm8;
if(isfield(mp.dm8,'V')==false); mp.dm8.V = mp.dm8.V0coef*ones(mp.dm8.NactTotal,1); else; mp.dm8.V = mp.DM8V0; end %--Initial DM voltages
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
mp.dm8.inf_datacube =  padOrCropEven(1-falco_gen_annular_FPM(FPMgenInputs),mp.dm8.NdmPad);
%--Make or read in DM8 disk for the compact model
FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
mp.dm8.compact.inf_datacube =  padOrCropEven(1-falco_gen_annular_FPM(FPMgenInputs),mp.dm8.compact.NdmPad);
%         figure(204); imagesc(mp.dm8.inf_datacube); axis xy equal tight; colorbar;
%         figure(205); imagesc(mp.dm8.compact.inf_datacube); axis xy equal tight; colorbar;



%%

%% OLD CODE
% 
% if(isfield(mp,'flagDM9inf3x3'))
%     mp.dm9.xcent_dm = mp.dm9.Nact/2 - 1/2;
%     mp.dm9.ycent_dm = mp.dm9.Nact/2 - 1/2;
%     if(strcmpi(mp.centering,'interpixel'))
%        error('falco_init_ws: The 3x3 influence function for DM9 requires a pixel-centered coordinate system (for now).');
%     end
% else
%     mp.dm9.xcent_dm = mp.dm9.Nact/2 - 1/2;
%     mp.dm9.ycent_dm = mp.dm9.Nact/2 - 1/2;
% end
% mp.dm9.dm_spacing = 1/mp.dm9.actres*(mp.fl*mp.lambda0/mp.P2.D); % meters, pitch of DM actuators
% mp.dm9.compact = mp.dm9;
% 
% mp.dm9.compact.dm_spacing = mp.dm9.dm_spacing; % meters, pitch of DM actuators
% 
% mp.dm9.dx_inf0 = (mp.dm9.dx_inf0_act)*mp.dm9.dm_spacing;%1e-4; % meters, sampling of the influence function 
% 
% mp.dm9.compact.dx_inf0 = (mp.dm9.compact.dx_inf0_act)*mp.dm9.compact.dm_spacing;%1e-4; % meters, sampling of the influence function 
% 
% mp.dm9.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.res; % width of a pixel at the FPM in the full model (meters)
% mp.dm9.compact.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.res; % width of a pixel at the FPM in the compact model (meters)
% 
% if(length(mp.dm9.inf0)==3)
%     mp.dm9 = falco_fpm_inf_cube_3x3(mp.dm9);
%     mp.dm9.compact = falco_fpm_inf_cube_3x3(mp.dm9.compact);
% else
%     mp.dm9 = falco_fpm_inf_cube(mp.dm9);
%     mp.dm9.compact = falco_fpm_inf_cube(mp.dm9.compact);
% end
%         
% %% Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)
% r_cent_lam0D = mp.dm9.r_cent_act*mp.dm9.dm_spacing/(mp.dm9.dxi)/mp.F3.full.res;
% 
% 
%%
% mp.F3.RinA_inds = [];
% mp.F3.RinAB_inds = [];
% for ii=1:mp.dm9.NactTotal
%     %--Zero out FPM actuators beyond the allowed radius (mp.F3.Rin)
%     if(r_cent_lam0D(ii) > mp.F3.Rin-mp.dm9.FPMbuffer)
%        mp.dm9.inf_datacube(:,:,ii) = 0*mp.dm9.inf_datacube(:,:,ii);
%        mp.dm9.compact.inf_datacube(:,:,ii) = zeros(size(mp.dm9.compact.inf_datacube(:,:,ii)));
%     end
% 
%     %--Get the indices for the actuators within radius mp.F3.RinA
%     if(r_cent_lam0D(ii) <= mp.F3.RinA-mp.dm9.FPMbuffer)
%         mp.F3.RinA_inds = [mp.F3.RinA_inds; ii];
%     else %--Get the indices for the actuators between radii mp.F3.RinA and mp.F3.Rin
%         mp.F3.RinAB_inds = [mp.F3.RinAB_inds; ii];
%     end
% end
% fprintf('%d actuators in DM9.\n',mp.dm9.NactTotal);
% 
% mp.dm9.ABfac = 1;%1/5; %--Gain factor between inner and outer FPM regions
% mp.dm9.VtoHavg = 1e-9; %--gain of DM9 (meters/Volt)
% mp.dm9.VtoH = mp.dm9.VtoHavg*ones(mp.dm9.NactTotal,1);%1*1e-9*ones(mp.dm9.Nact); % Gains: volts to meters in surface height;
% mp.dm9.VtoH(mp.F3.RinAB_inds) = mp.dm9.ABfac*mp.dm9.VtoH(mp.F3.RinAB_inds);
% 
% if(isfield(mp.dm9,'V')==false); mp.dm9.V = zeros(mp.dm9.NactTotal,1); mp.dm9.V(mp.F3.RinA_inds) = mp.dm9.V0coef; else; mp.dm9.V = mp.DM9V0; end %--Initial DM9 voltages
% 
% mp.dm9.Vmin = min(mp.t_diel_nm_vec);  % minimum thickness of FPM dielectric layer (nm)
% mp.dm9.Vmax = max(mp.t_diel_nm_vec);  % maximum thickness (from one actuator, not of the facesheet) of FPM dielectric layer (nm)

%%
% %%%%%---OPTIONS FOR DEFINING DM8 (FPM Metal)
% mp.dm8.VtoHavg = 1e-9; %--gain of DM8 (meters/Volt)
% mp.dm8.Vmin = min(mp.t_metal_nm_vec); % minimum thickness of FPM metal layer (nm)
% mp.dm8.Vmax = max(mp.t_metal_nm_vec); % maximum thickness (from one actuator, not of the facesheet) of FPM metal layer (nm)
% 
% %         %--DM8 Option 1: Set as same basis set (DM actuators) as DM9. 
% %         dm8v0coef = mp.dm8.V0coef; %--Starting value. Save as temporary value to avoid overwriting
% %         mp.dm8 = mp.dm9; %--DANGER: OVERWRITING ANY PREVIOUS mp.dm8 DATA!!!
% %         mp.dm8.V0coef = dm8v0coef;
% %         fprintf('%d actuators in DM8.\n',mp.dm8.NactTotal);
% %         mp.dm8.ABfac = 1;%1/5; %--Gain factor between inner and outer FPM regions
% %         mp.dm8.VtoH = mp.dm8.VtoHavg*ones(mp.dm8.NactTotal,1);%1*1e-9*ones(mp.dm9.Nact); % Gains: volts to meters in surface height;
% %         mp.dm8.VtoH(mp.F3.RinAB_inds) = mp.dm8.ABfac*mp.dm8.VtoH(mp.F3.RinAB_inds);
% %         if(isfield(mp.dm8,'V')==false); mp.dm8.V = zeros(mp.dm8.NactTotal,1); mp.dm8.V(mp.F3.RinA_inds) = mp.dm8.V0coef; else; mp.dm8.V = DM8V0; end %--Initial DM voltages
% 
% 
% %--DM8 Option 2: Set basis as a single nickel disk.
% mp.dm8.NactTotal = 1;
% mp.dm8.act_ele = 1;
% fprintf('%d actuators in DM8.\n',mp.dm8.NactTotal);
% mp.dm8.VtoH = mp.dm8.VtoHavg*ones(mp.dm8.NactTotal,1); % Gains: volts to meters in surface height;
% mp.dm8.xy_box_lowerLeft = [1;1];
% mp.dm8.compact = mp.dm8;
% if(isfield(mp.dm8,'V')==false); mp.dm8.V = mp.dm8.V0coef*ones(mp.dm8.NactTotal,1); else; mp.dm8.V = mp.DM8V0; end %--Initial DM voltages
% % Don't define extra actuators and time:
% if(mp.F3.Rin~=mp.F3.RinA); error('falco_init_ws.m: Change mp.F3.Rin and mp.F3.RinA to be equal to avoid wasting time.'); end
% 
% % Copy over some common values from DM9:
% mp.dm8.dxi = mp.dm9.dxi; %--Width of a pixel at the FPM in full model (meters)
% mp.dm8.NdmPad = mp.dm9.NdmPad;
% mp.dm8.Nbox = mp.dm8.NdmPad;
% mp.dm8.compact.dxi = mp.dm9.compact.dxi; %--Width of a pixel at the FPM in compact model (meters)
% mp.dm8.compact.NdmPad = mp.dm9.compact.NdmPad;
% mp.dm8.compact.Nbox = mp.dm8.compact.NdmPad;
% 
% %--Make or read in DM8 disk for the full model
% FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
% FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
% FPMgenInputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
% FPMgenInputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
% FPMgenInputs.centering = mp.centering;
% mp.dm8.inf_datacube =  padOrCropEven(1-falco_gen_annular_FPM(FPMgenInputs),mp.dm8.NdmPad);
% %--Make or read in DM8 disk for the compact model
% FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
% mp.dm8.compact.inf_datacube =  padOrCropEven(1-falco_gen_annular_FPM(FPMgenInputs),mp.dm8.compact.NdmPad);
% %         figure(204); imagesc(mp.dm8.inf_datacube); axis xy equal tight; colorbar;
% %         figure(205); imagesc(mp.dm8.compact.inf_datacube); axis xy equal tight; colorbar;
% 
% %--Zero out parts of DM9 actuators that go outside the nickel disk. Also apply the grayscale edge.
% DM8windowFull = mp.dm8.inf_datacube;
% DM8windowCompact = mp.dm8.compact.inf_datacube;
% %         DM8windowFull = mp.dm8.inf_datacube > 0;
% %         DM8windowCompact = mp.dm8.compact.inf_datacube > 0;
% for iact=1:mp.dm9.NactTotal 
%     if( any(any(mp.dm9.inf_datacube(:,:,iact))) && any(mp.dm9.VtoH(iact)) )
%         %fprintf('iact = %d\n',iact);
%         y_box_ind = mp.dm9.xy_box_lowerLeft(1,iact):mp.dm9.xy_box_lowerLeft(1,iact)+mp.dm9.Nbox-1; % x-indices in pupil arrays for the box
%         x_box_ind = mp.dm9.xy_box_lowerLeft(2,iact):mp.dm9.xy_box_lowerLeft(2,iact)+mp.dm9.Nbox-1; % y-indices in pupil arrays for the box
%         mp.dm9.inf_datacube(:,:,iact) = DM8windowFull(x_box_ind,y_box_ind).*mp.dm9.inf_datacube(:,:,iact);
% 
%         y_box_ind = mp.dm9.compact.xy_box_lowerLeft(1,iact):mp.dm9.compact.xy_box_lowerLeft(1,iact)+mp.dm9.compact.Nbox-1; % x-indices in pupil arrays for the box
%         x_box_ind = mp.dm9.compact.xy_box_lowerLeft(2,iact):mp.dm9.compact.xy_box_lowerLeft(2,iact)+mp.dm9.compact.Nbox-1; % y-indices in pupil arrays for the box
%         mp.dm9.compact.inf_datacube(:,:,iact) = DM8windowCompact(x_box_ind,y_box_ind).*mp.dm9.compact.inf_datacube(:,:,iact);
%     end
% end




end %--END OF FUNCTION
