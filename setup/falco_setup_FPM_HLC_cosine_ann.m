% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% REVISION HISTORY:
% --------------
% Modified on 2019-06-13 by A.J. Riggs from falco_setup_FPM_HLC_3foldZern.m
% to falco_setup_FPM_HLC_cosine.m
% Created by A.J. Riggs on 2018-10-01 by extracting material from
% falco_init_ws.m.
% ---------------

function mp = falco_setup_FPM_HLC_cosine(mp)

%% DM8 and DM9 (Optimizable FPM) Setup

%--Centering of DM surfaces on array
mp.dm8.centering = mp.centering;
mp.dm9.centering = mp.centering;

mp.dm9.compact = mp.dm9;

%% DM9 as Cosine Rings

% <<< DEBUGGING: HARD-CODED VALUES FOR TESTING
mp.dm9.actres = 5;
mp.F3.Rin = 2.8;
mp.F3.compact.res = 100;
mp.F3.full.res = 50;
mp.centering = 'pixel';
mp.dm9.VtoHavg = 1e-9;
mp.fl = 1;
mp.lambda0 = 1;
mp.P2.D = 1;
% >>>DEBUGGING


%--Pixel size [meters]
mp.dm9.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.res; % width of a pixel at the FPM in the full model (meters)
mp.dm9.compact.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.res; % width of a pixel at the FPM in the compact model (meters)


drCos = 1/mp.dm9.actres; %--Width and double-separation of the cosine rings [lambda0/D]
mp.dm9.NactTotal = ceil(2*mp.dm9.actres*mp.F3.Rin);

%--Generate datacube of influence functions, which are rings with radial cosine profile
%--Compact model
NbeamCompact = mp.F3.compact.res*mp.F3.Rin*2;
% mp.dm9.compact.inf_datacube = falco_3fold_symmetry_Zernikes(NbeamCompact,mp.dm9.maxRadialOrder,mp.centering,'SymmAxis','y');
mp.dm9.compact.NdmPad = ceil_even(1+2*mp.F3.Rin*mp.F3.compact.res);% size(mp.dm9.compact.inf_datacube,1);
mp.dm9.compact.Nbox = mp.dm9.compact.NdmPad; %--the modes take up the full array.
%--Normalized coordinates: Compact model
if(strcmpi(mp.centering,'pixel')  ) 
    xc = (-mp.dm9.compact.NdmPad/2:(mp.dm9.compact.NdmPad/2-1))/mp.F3.compact.res;    
else
    xc = (-(mp.dm9.compact.NdmPad-1)/2:(mp.dm9.compact.NdmPad-1)/2-1)/mp.F3.compact.res;    
end
[Xc,Yc] = meshgrid(xc);
Rc = sqrt(Xc.^2 + Yc.^2);
mp.dm9.compact.inf_datacube = zeros(mp.dm9.compact.NdmPad,mp.dm9.compact.NdmPad,mp.dm9.NactTotal);

hg_expon = 44; %--Found empirically
apRad = mp.F3.Rin/(mp.F3.Rin+0.1); %--Found empirically
OD = 1;
mask = Rc<=mp.F3.Rin;
windowFull = mask.*exp(-(Rc/mp.F3.Rin/(apRad*OD)).^hg_expon);

drSep = drCos/2;
%--Compute the ring influence functions
for ri = 1:mp.dm9.NactTotal
    modeTemp = windowFull.*(1+(-1)^(mod(ri+1,2))*cos(2*pi*(Rc*mp.dm9.actres-0.5)))/2;
    rMin = drSep*(ri - 1);
    rMax = drSep*(ri + 1);
    if(ri==1) %--Fill in the center
        modeTemp(Rc<drSep) = 1; 
    else
        modeTemp(Rc<rMin) = 0;
    end
    modeTemp(Rc>rMax) = 0;
    mp.dm9.compact.inf_datacube(:,:,ri) = modeTemp;
    figure(10); imagesc(xc,xc,mp.dm9.compact.inf_datacube(:,:,ri)); axis xy equal tight; colorbar; drawnow;
    figure(11); plot(xc,mp.dm9.compact.inf_datacube(:,mp.dm9.compact.NdmPad/2+1,ri)); xlim([0,mp.F3.Rin]); drawnow;
    pause(0.1);
end
figure(12); imagesc(xc,xc,sum(mp.dm9.compact.inf_datacube,3)); axis xy equal tight; colorbar; drawnow;



% for ri = 1:mp.dm9.NactTotal/2   
%     Rcenter = drCos*(ri - 0.25);
%     modeTemp = windowFull.*(1+cos(2*pi.*(Rc*mp.dm9.actres-1/2)));
%     rMin = drCos*(ri - 1.);
%     rMax = drCos*(ri - 0);
%     modeTemp(Rc<rMin) = 0;
%     modeTemp(Rc>rMax) = 0;
%     mp.dm9.compact.inf_datacube(:,:,ri) = modeTemp;
%     figure(10); imagesc(xc,xc,mp.dm9.compact.inf_datacube(:,:,ri)); axis xy equal tight; colorbar; drawnow;
% %     figure(11); plot(xc,mp.dm9.compact.inf_datacube(:,mp.dm9.compact.NdmPad/2+1,ri)); xlim([0,mp.F3.Rin]); drawnow;
%     pause(0.1);
% end
% 
%     figure(12); imagesc(xc,xc,sum(mp.dm9.compact.inf_datacube,3)); axis xy equal tight; colorbar; drawnow;


%%


% %--Full Model
% NbeamFull = mp.F3.full.res*mp.F3.Rin*2;
% % mp.dm9.inf_datacube = falco_3fold_symmetry_Zernikes(NbeamFull,mp.dm9.maxRadialOrder,mp.centering,'SymmAxis','y');
% mp.dm9.NdmPad = ceil_even(1+2*mp.F3.Rin*mp.F3.res);%size(mp.dm9.inf_datacube,1);
% mp.dm9.Nbox = mp.dm9.NdmPad; %--the modes take up the full array.
% %--Normalized coordinates: Full model
% if(strcmpi(mp.centering,'pixel')  ) 
%     xf = (-mp.dm9.NdmPad/2:(mp.dm9.NdmPad/2-1))/mp.F3.res;    
% else
%     xf = (-(mp.dm9.NdmPad-1)/2:(mp.dm9.NdmPad-1)/2-1)/mp.F3.res;    
% end
% [Xf,Yf] = meshgrid(xf);
% Rf = sqrt(Xf.^2 + Yf.^2);
% 
% 
% mp.dm9.inf_datacube = zeros(mp.dm9.NdmPad,mp.dm9.NdmPad,mp.dm9.NactTotal);

mp.dm9.NactTotal = size(mp.dm9.inf_datacube,3);
mp.dm9.VtoH = mp.dm9.VtoHavg*ones(mp.dm9.NactTotal,1);

%--Lower-left pixel coordinates are all (1,1) since the Zernikes take up the full array.
mp.dm9.xy_box_lowerLeft = ones(2,mp.dm9.NactTotal);
mp.dm9.compact.xy_box_lowerLeft = ones(2,mp.dm9.NactTotal);

%--Coordinates for the full FPM array
if(strcmpi(mp.centering,'pixel')  ) 
    mp.dm9.compact.x_pupPad = (-mp.dm9.compact.NdmPad/2:(mp.dm9.compact.NdmPad/2-1))*mp.dm9.compact.dxi; % meters, coords for the full DM arrays. Origin is centered on a pixel
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

mp.dm9.Vmin = min(mp.t_diel_nm_vec); % minimum thickness of FPM dielectric layer (nm)
mp.dm9.Vmax = max(mp.t_diel_nm_vec); % maximum thickness (from one actuator, not of the facesheet) of FPM dielectric layer (nm)

%% Apply window to DM9 influence functions
% 
% hg_expon = 44; %--Found empirically
% apRad = mp.F3.Rin/(mp.F3.Rin+0.1); %--Found empirically
% OD = 1;
% 
% %--Full Model
% Narray = mp.dm9.NdmPad;
% Nbeam = NbeamFull;
% %--Coordinates normalized to the beam radius (not diameter)
% switch mp.centering
%     case 'interpixel'
%         xs = (-(Narray-1)/2:(Narray-1)/2)/Nbeam*2;
%     case 'pixel'
%         xs = (-Narray/2:(Narray/2-1))/Nbeam*2;
% end
% [XS,YS] = meshgrid(xs);
% RS = sqrt(XS.^2+YS.^2); 
% mask = RS<=1;
% windowFull = mask.*exp(-(RS/(apRad*OD)).^hg_expon);
% 
% %--Compact Model
% Narray = mp.dm9.compact.NdmPad;
% Nbeam = NbeamCompact;
% %--Coordinates normalized to the beam radius (not diameter)
% switch mp.centering
%     case 'interpixel'
%         xs = (-(Narray-1)/2:(Narray-1)/2)/Nbeam*2;
%     case 'pixel'
%         xs = (-Narray/2:(Narray/2-1))/Nbeam*2;
% end
% [XS,YS] = meshgrid(xs);
% RS = sqrt(XS.^2+YS.^2); 
% mask = RS<=1;
% windowCompact = mask.*exp(-(RS/(apRad*OD)).^hg_expon);
% 
% 
% for di = 1:mp.dm9.NactTotal
%     mp.dm9.inf_datacube(:,:,di) = windowFull.*mp.dm9.inf_datacube(:,:,di);
%     mp.dm9.compact.inf_datacube(:,:,di) = windowCompact.*mp.dm9.compact.inf_datacube(:,:,di);
% end

%% DM8 as a disk (piston as only influence function)

%%%%%---OPTIONS FOR DEFINING DM8 (FPM Metal)
mp.dm8.VtoHavg = 1e-9; %--gain of DM8 (meters/Volt) %--MOVE TO overwritable defaults when controlling DM8.
mp.dm8.Vmin = min(mp.t_metal_nm_vec); % minimum thickness of FPM metal layer (nm)
mp.dm8.Vmax = max(mp.t_metal_nm_vec); % maximum thickness (from one actuator, not of the facesheet) of FPM metal layer (nm)

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
mp.dm8.inf_datacube = padOrCropEven(1-falco_gen_annular_FPM(FPMgenInputs),mp.dm8.NdmPad);
%--Make or read in DM8 disk for the compact model
FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
mp.dm8.compact.inf_datacube = padOrCropEven(1-falco_gen_annular_FPM(FPMgenInputs),mp.dm8.compact.NdmPad);

%%

end %--END OF FUNCTION