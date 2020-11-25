% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Function to generate Roman CGI pupil 20200513 in Matlab using PROPER
%  functions for rectangles and ellipses.
%
%  Parameters were found originally found by Dwight Moody and then adjusted
%  by A.J. Riggs by fitting ellipses and rectangles to the universally
%  released file 'pupil_CGI-20200513_8k_binary_noF.png'.
%
% %--Order of operations for pupil modifications:
% 1) Pad
% 2) Magnify
% 3) Rotate
% 4) Translate
%
% NOTE: All pupil features have normalized length units of pupil diameters.

function pupil = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, varargin)

%% Define the best-fit values for ellipses and rectangles (DO NOT CHANGE)

primaryRadiusYpixels = 4027.25;
ODpixels = 2.*primaryRadiusYpixels;

primaryRadiusX = 3990.0/ODpixels;
primaryRadiusY = primaryRadiusYpixels/ODpixels;
primaryCenterX = 0.;
primaryCenterY = 0.;

secondaryRadiusX = 1209.65/ODpixels;
secondaryRadiusY = 1220.0/ODpixels;
secondaryCenterX = 0.0/ODpixels;
secondaryCenterY = -2.95/ODpixels;

strutEndVecX1 = ([843.9, 728.0, 47.5, -192.0, -676.0, -816.65])/ODpixels;
strutEndVecY1 = ([550.85, 580.35, -970.65, -1097.15, 605.55, 458.85])/ODpixels;

strutEndVecX2 = ([1579.9, 3988.0, 2429.7, -2484.0, -3988.0, -1572.65])/ODpixels;
strutEndVecY2 = ([3866.85, -511.65, -3214.65, -3256.15, -504.45, 3866.85])/ODpixels;

strutCenterVecX = (strutEndVecX1 + strutEndVecX2)/2.0;
strutCenterVecY = (strutEndVecY1 + strutEndVecY2)/2.0;

strutWidthVec = [257.0, 259.0, 258.0, 258.0, 259.0, 257.0]/ODpixels;

strutAngleVec = atan2(strutEndVecY2-strutEndVecY1, strutEndVecX2-strutEndVecX1)*(180/pi); % [degrees]

tabRadiusVecX = [1343.0, 1343.0, 1364.0]/ODpixels;
tabRadiusVecY = [1353.0, 1353.0, 1374.0]/ODpixels;
tabCenterVecX = [0.0, 0.0, 0.0]/ODpixels;
tabCenterVecY = [53.85, 53.85, 67.6]/ODpixels;

lStrut = 0.55;

deltaAngle = 2.5*pi/16;
angTabStart = [0.616 - deltaAngle/2.0; 2.54 - deltaAngle/2.0; -1.57 - deltaAngle/2.0];
angTabEnd   = [0.616 + deltaAngle/2.0; 2.54 + deltaAngle/2.0; -1.57 + deltaAngle/2.0];

xcStrutVec = strutCenterVecX;
ycStrutVec = strutCenterVecY;
angStrutVec = strutAngleVec;
wStrutVec = strutWidthVec;

%% Changes to the pupil

if(size(varargin,2)>1)
    error('falco_gen_pupil_Roman_CGI_20200513.m: Too many inputs')
elseif(size(varargin,2)==1)
    changes = varargin{1}; %--Structure containing which values to change;
else
    changes.dummy = 1; %--Else initialize structure
end

%%--(Optional) Lyot stop mode (concentric, circular ID and OD)
if(~isfield(changes,'flagLyot'))
    flagLyot = false; 
else
    flagLyot = changes.flagLyot; 
end
if(flagLyot == true)
    if(isfield(changes,'ID'))
        ID = changes.ID;
    else
        error('changes.ID must be defined for Lyot stop generation mode.')
    end

    if(isfield(changes,'OD'))
        OD = changes.OD;
    else
        error('changes.OD must be defined for Lyot stop generation mode.')
    end
end

%%--Oversized strut features: overwrite defaults if values specified.
if(isfield(changes,'wStrut')); wStrutVec = changes.wStrut*ones(6,1);  end
if(isfield(changes,'wStrutVec')); wStrutVec = changes.wStrutVec;  end %--Unlikely to be called. Would overrule line above this one.

%%--Padding values for obscurations
%--Defaults of Bulk Changes: (All length units are pupil diameters. All angles are in degrees.)
if(~isfield(changes,'xShear')); changes.xShear = 0.; end
if(~isfield(changes,'yShear')); changes.yShear = 0.; end
if(~isfield(changes,'magFac')); changes.magFac = 1.0; end
if(~isfield(changes,'clock_deg')); changes.clock_deg = 0.0; end
if(~isfield(changes,'flagRot180')); changes.flagRot180 = false; end

%--Defaults for obscuration padding: (All length units are pupil diameters.)
if(~isfield(changes,'pad_all')); changes.pad_all = 0.0; end
if(~isfield(changes,'pad_strut')); changes.pad_strut = 0.0; end
if(~isfield(changes,'pad_COBS')); changes.pad_COBS = 0.0; end
if(~isfield(changes,'pad_COBStabs')); changes.pad_COBStabs = 0.0; end
if(~isfield(changes,'pad_OD')); changes.pad_OD = 0.0; end

%--Values to use for bulk clocking, magnification, and translation
xShear = changes.xShear;
yShear = changes.yShear;
magFac = changes.magFac;
clock_deg = changes.clock_deg;
clock_rad = deg2rad(clock_deg);
flagRot180 = changes.flagRot180;

%--Padding values. (pad_all is added to all the rest)
pad_all = changes.pad_all; %--Uniform padding on all features
pad_strut = changes.pad_strut + pad_all; % amount padded to each side of strut
pad_COBS = changes.pad_COBS + pad_all;
pad_COBStabs = changes.pad_COBStabs + pad_all;
pad_OD = changes.pad_OD + pad_all; %--Radial padding at the edge

%--Rotation matrix used on center coordinates.
rotMat = [cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)];

%% Nominal Mask Coordinates

%--number of points across output array.
if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam*max([magFac, 1.0]) + 1 + 2*magFac*Nbeam*max(abs([xShear, yShear]))); % Sometimes requires two more pixels when pixel centered.
elseif(strcmpi(centering,'interpixel'))
    Narray = ceil_even(Nbeam*max([magFac, 1.0]) + 2*magFac*Nbeam*max(abs([xShear, yShear]))); %--No zero-padding needed if beam is centered between pixels
end
if(isfield(changes,'Narray')); Narray = changes.Narray; end

if(strcmpi(centering,'interpixel'))
    xs = (-(Narray-1)/2:(Narray-1)/2)/Nbeam;
elseif(strcmpi(centering,'pixel'))
    xs = (-(Narray/2):Narray/2-1)/Nbeam;
end

[XS,YS] = meshgrid(xs);

%% PROPER SETUP VALUES

Dbeam = 1;                 %--Diameter of aperture, normalized to itself
wl   = 1e-6;              % wavelength (m); Dummy value--no propagation here, so not used.
bdf = Nbeam/Narray; %--beam diameter factor in output array
dx = Dbeam/Nbeam;

switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; % = -dx/2/bdf;
    case {'pixel'}
        cshift = 0;
        if(flagRot180)
            cshift = -dx; % = -dx/bdf;
        end
end

%% INITIALIZE PROPER Wave Structure for Struts
bm = prop_begin(Dbeam, wl, Narray,'beam_diam_fraction',bdf);
nSubsamples = 101;
prop_set_antialiasing(nSubsamples)

%--Struts
for iStrut=1:6
    angDeg = angStrutVec(iStrut) + clock_deg; % [degrees]
    wStrut = magFac*(wStrutVec(iStrut) + 2*pad_strut);
    lStrutIn = magFac*lStrut;
    xc = magFac*(xcStrutVec(iStrut));
    yc = magFac*(ycStrutVec(iStrut));
    cxy = rotMat*[xc; yc];
    xc = cxy(1)+xShear;
    yc = cxy(2)+yShear;
    bm = prop_rectangular_obscuration(bm, lStrutIn, wStrut, 'XC',xc+cshift, 'YC',yc+cshift, 'ROTATION',angDeg);
end

if(flagLyot == false)
    %% PRIMARY MIRROR (OUTER DIAMETER)
    ra_OD_x = magFac*(primaryRadiusX-pad_OD);
    ra_OD_y = magFac*(primaryRadiusY-pad_OD);
    cx_OD = magFac*primaryCenterX;
    cy_OD = magFac*primaryCenterY;
    cxy = rotMat*[cx_OD; cy_OD];
    cx_OD = cxy(1) + xShear + cshift;
    cy_OD = cxy(2) + yShear + cshift;
    bm = prop_elliptical_aperture( bm, ra_OD_x, ra_OD_y, 'XC', cx_OD, 'YC', cy_OD, 'ROTATION', -clock_deg);

    %% SECONDARY MIRROR (INNER DIAMETER)
    ra_ID_x = magFac*(secondaryRadiusX + pad_COBS);
    ra_ID_y = magFac*(secondaryRadiusY + pad_COBS);
    cx_ID = magFac*secondaryCenterX;
    cy_ID = magFac*secondaryCenterY;
    cxy = rotMat*[cx_ID; cy_ID];
    cx_ID = cxy(1) + xShear + cshift;
    cy_ID = cxy(2) + yShear + cshift;
    bm = prop_elliptical_obscuration(bm, ra_ID_x, ra_ID_y,'XC',cx_ID,'YC',cy_ID,'ROTATION',-clock_deg);

    %--Tabs where Struts Meet Central Obscuration (aka COBS)
    nTabs = 3;
    tabCube = ones(Narray, Narray, nTabs);
    
    for iTab = 1:nTabs
        cobsTabsMask = zeros(Narray);
        
        xyShear = rotMat*[tabCenterVecX(iTab); tabCenterVecY(iTab)];

        XSnew = XS + magFac*xyShear(1) - xShear;
        YSnew = YS + magFac*xyShear(2) - yShear;
        THETAS = atan2(YSnew, XSnew);

        ang1 = angTabStart(iTab) + clock_rad;
        ang2 = angTabEnd(iTab) + clock_rad;
        while (ang2 > 2*pi) && (ang1 > 2*pi)
            ang1 = ang1 - 2*pi;
            ang2 = ang2 - 2*pi;
        end
        while (ang2 < -2*pi) && (ang1 < -2*pi)
            ang1 = ang1 + 2*pi;
            ang2 = ang2 + 2*pi;
        end
        if ((ang2 < 2*pi) && (ang2 > pi)) && ((ang1 < pi) || ang1 > 0)
            THETAS(THETAS<0) = THETAS(THETAS<0) + 2*pi;   
        elseif (ang2 > 2*pi) && ((ang1 > pi) && (ang1 < 2*pi))
            THETAS = THETAS + 2*pi; 
        elseif (ang1 < -pi) && (ang2 > -pi)
            THETAS(THETAS>0) = THETAS(THETAS>0) - 2*pi;
        elseif (ang1 < -pi) && (ang2 < -pi)
            THETAS = THETAS - 2*pi; 
        end
        cobsTabsMask( (THETAS <= ang2) & (THETAS >= ang1)) = 1.0;

        %--ELLIPSE:

        % Full ellipse to be multiplied by the mask to get just tabs
        bm2 = prop_begin(Dbeam, wl, Narray,'beam_diam_fraction',bdf);
        prop_set_antialiasing(nSubsamples)
        cx_tab = magFac*tabCenterVecX(iTab);
        cy_tab = magFac*tabCenterVecY(iTab);
        cxy = rotMat*[cx_tab; cy_tab];
        cx_tab = cxy(1)+xShear;
        cy_tab = cxy(2)+yShear;
        tabRadiusX = magFac*(tabRadiusVecX(iTab) + pad_COBStabs);
        tabRadiusY = magFac*(tabRadiusVecY(iTab) + pad_COBStabs);
        bm2 = prop_elliptical_obscuration(bm2, tabRadiusX, tabRadiusY,'XC',cx_tab+cshift,'YC',cy_tab+cshift,'ROTATION',-clock_deg);
        tabEllipse = 1-ifftshift(abs(bm2.wf));

        tabSector = cobsTabsMask.*tabEllipse;

        tabCube(:,:,iTab) = 1-tabSector;
%         figure(11); imagesc(cobsTabsMask); axis xy equal tight; colorbar; drawnow;
%         figure(12); imagesc(tabEllipse); axis xy equal tight; colorbar; drawnow;
%         figure(13); imagesc(tabCube(:,:,iTab)); axis xy equal tight; colorbar; drawnow;
    end

    pupil = tabCube(:,:,1).*tabCube(:,:,2).*tabCube(:,:,3).*ifftshift(abs(bm.wf));

else % (Lyot stop mode)

    %--OUTER DIAMETER
    ra_OD = magFac*(OD/2.);
    cx_OD = xShear;
    cy_OD = yShear;
    bm = prop_circular_aperture(bm, ra_OD, 'XC', cx_OD+cshift, 'YC', cy_OD+cshift);

    %--INNER DIAMETER
    ra_ID = magFac*(ID/2.);
    cx_ID = xShear;
    cy_ID = yShear;
    bm = prop_circular_obscuration(bm, ra_ID, 'XC', cx_ID+cshift, 'YC', cy_ID+cshift);

    pupil = ifftshift(abs(bm.wf));

end

if(flagRot180)
   pupil = rot90(pupil,2);
end

end % EOF
