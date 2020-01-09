% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Function to generate WFIRST pupil 20191009 in Matlab using PROPER 
%    functions for rectangles and FALCO for ellipses.
%
%  Parameters were found by Dwight Moody by fitting ellipses and rectangles to the 
%    universally released file '2019-10-09 CGI entrance pupil 8k binary.png'.
%
% %--Order of Operations:
% 1) Pad
% 2) Magnify
% 3) Rotate
% 4) Translate
%
% NOTE: All pupil features have normalized length units of pupil diameters.

function pupil = falco_gen_pupil_WFIRST_CGI_20191009(Nbeam,centering,varargin)

%% Define the best-fit values for ellipses and rectangles (DO NOT CHANGE)

primaryRadiusYpixels = 4022.0;
ODpixels = 2*primaryRadiusYpixels;

primaryRadiusX = 3987.0/ODpixels;
primaryRadiusY = primaryRadiusYpixels/ODpixels;
primaryCenterX = 0.0;
primaryCenterY = 0.0;
% primarradiusYYpixels = 0*primarradiusYY*ODpixels;

secondarradiusYX = 1208.5/ODpixels;
secondarradiusYY = 1219.0/ODpixels;
% secondaryCenterX = -0.5/ODpixels;
% secondaryCenterY = -1.3/ODpixels;
secondaryCenterX = -0.0/ODpixels;
secondaryCenterY = -1.8/ODpixels;

strutEndVecX1 = ([843.0, 728.0, 47.0, -192.0, -676.0, -816.0])/ODpixels;
strutEndVecY1 = ([552.0, 581.0, -968.0, -1095.0, 606.5, 460.0])/ODpixels;

strutEndVecX2 = ([1579.0, 3988.0, 2430.0, -2484.0, -3988.0, -1572.0])/ODpixels;
strutEndVecY2 = ([3868.0, -511.0, -3212.0, -3254.0, -503.5, 3868.0])/ODpixels;

strutCenterVecX = (strutEndVecX1 + strutEndVecX2)/2.0;
strutCenterVecY = (strutEndVecY1 + strutEndVecY2)/2.0;

strutWidthVec = [257.0, 259.0, 258.0, 258.0, 258.5, 257.0]/ODpixels;

strutAngleVec = atan2(strutEndVecY2-strutEndVecY1, strutEndVecX2-strutEndVecX1)*(180/pi);

tabRadiusVecX = [1342.0, 1342.0, 1364.0]/ODpixels;
tabRadiusVecY = [1352.0, 1352.0, 1374.0]/ODpixels;
tabCenterVecX = [0.0, 0.0, 0.0]/ODpixels;
% tabCenterVecY = [60.0, 60.0, 60.0]/ODpixels;
tabCenterVecY = [55.0, 55.0, 70.0]/ODpixels;

lStrut = 0.55;

deltaAngle = 2.5*pi/16;
angTabStart = [0.616 - deltaAngle/2.0; 2.54 - deltaAngle/2.0; -1.57 - deltaAngle/2.0]; 
angTabEnd   = [0.616 + deltaAngle/2.0; 2.54 + deltaAngle/2.0; -1.57 + deltaAngle/2.0]; 

xcStrutVec = strutCenterVecX;
ycStrutVec = strutCenterVecY;
angStrutVec = strutAngleVec;
wStrutVec = strutWidthVec;

%  #### primary #####
%   dyn1.add_circ(ry=4022.5, rx=3987.0, shifty=0, circ_invert=False, aper_invert=False, tag='a1')
% 
%   #### secondary #####
%   dyn1.add_circ(ry=1219., rx=1208.5, shifty=-1.3, shiftx=-0.5, circ_invert=True, aper_invert=False, tag='a2')
% 
%   ######################################### struts #####
%   dyn1.add_blockline(tag="a3",y1=552.0,x1=843.0,y2=3868.0,x2=1579.0,width=257.0, \
%       block_type=1,aper_join_type=2,insert_i=-1,insert_after=True)
%   ####
%   dyn1.add_blockline(tag="a4",y1=581.0,x1=728.0,y2=-511.0,x2=3988.0,width=259.0, \
%       block_type=1,aper_join_type=2,insert_i=-1,insert_after=True)
%   ####
%   dyn1.add_blockline(tag="a5",y1=-968.0,x1=47.0,y2=-3212.0,x2=2430.0,width=258.0, \
%       block_type=1,aper_join_type=2,insert_i=-1,insert_after=True)
%   ####
%   dyn1.add_blockline(tag="a6",y1=-1095.0,x1=-192.0,y2=-3254.0,x2=-2484.0,width=258.0, \
%       block_type=1,aper_join_type=2,insert_i=-1,insert_after=True)
%   ####
%   dyn1.add_blockline(tag="a7",y1=606.5,x1=-676.0,y2=-503.5,x2=-3988.0,width=258.5, \
%       block_type=1,aper_join_type=2,insert_i=-1,insert_after=True)
%   ####
%   dyn1.add_blockline(tag="a8",y1=460.0,x1=-816.0,y2=3868.0,x2=-1572.0,width=257.0, \
%       block_type=1,aper_join_type=2,insert_i=-1,insert_after=True)
% 
%   #### tabs #####
%   dyn1.add_circ(ry=1352, rx=1342, shifty=60.0, shiftx=0.0, angle=pi/6, block_angle=pi/16, tag='a9')
%   ###
%   dyn1.add_circ(ry=1352, rx=1342, shifty=60.0, shiftx=0.0, angle=2*pi/3+pi/6, block_angle=pi/16, tag='a10')
%   ###
%   dyn1.add_circ(ry=1374, rx=1364, shifty=60.0, shiftx=0.0, angle=4*pi/3+pi/6, block_angle=pi/16, tag='a11')


%% Changes to the pupil

if(size(varargin,2)>1)
    error('falco_gen_pupil_WFIRST_CGI_20191009.m: Too many inputs')
elseif(size(varargin,2)==1)
    changes = varargin{1}; %--Structure containing which values to change;
else
    changes.dummy = 1; %--Else initialize structure
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
xShear = changes.xShear;% - xcOD;
yShear = changes.yShear;% - ycOD;
magFac = changes.magFac;
clock_deg = changes.clock_deg;
clock_rad = deg2rad(clock_deg);
flagRot180 = changes.flagRot180;

%--Padding values. (pad_all is added to all the rest)
pad_all = changes.pad_all;%0.2/100; %--Uniform padding on all features
pad_strut = changes.pad_strut + pad_all;
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
else
    xs = (-(Narray/2):Narray/2-1)/Nbeam; 
end

[XS,YS] = meshgrid(xs);

%% PROPER SETUP VALUES

Dbeam = 1;                 %--Diameter of aperture, normalized to itself
wl   = 1e-6;              % wavelength (m); Dummy value--no propagation here, so not used.
bdf = Nbeam/Narray; %--beam diameter factor in output array
dx = Dbeam/Nbeam;

switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel','even'}
        cshift = -dx/2; % = -dx/2/bdf; 
    case {'pixel','odd'}
        cshift = 0;
        if(flagRot180)
            cshift = -dx; % = -dx/bdf; 
        end
end

%% INITIALIZE PROPER Wave Structure for Struts
bm = prop_begin(Dbeam, wl, Narray,'beam_diam_fraction',bdf);

%% PRIMARY MIRROR (OUTER DIAMETER)
ra_OD_x = magFac*(primaryRadiusX-pad_OD);
ra_OD_y = magFac*(primaryRadiusY-pad_OD);
cx_OD = magFac*primaryCenterX;
cy_OD = magFac*primaryCenterY;
cxy = rotMat*[cx_OD; cy_OD];
cx_OD = cxy(1) + xShear;% + cshift; --> cshift not needed because falco_gen_ellipse takes centering into account
cy_OD = cxy(2) + yShear;% + cshift; --> cshift not needed because falco_gen_ellipse takes centering into account
% bm = prop_elliptical_aperture( bm, ra_OD_x, ra_OD_y, 'XC', cx_OD, 'YC', cy_OD);

inputs.Narray = Narray;
inputs.Nbeam = Nbeam;
inputs.radiusX = ra_OD_x;
inputs.radiusY = ra_OD_y;
inputs.clockingDegrees = clock_deg;
inputs.centering = centering; 
inputs.xShear = cx_OD;
inputs.yShear = cy_OD;
% inputs.magFac = magFac; %  not needed because falco_gen_ellipse takes magnification into account
primaryAperture = falco_gen_ellipse(inputs);
% primaryAperture = 1;

%% SECONDARY MIRROR (INNER DIAMETER)
ra_ID_x = magFac*(secondarradiusYX + pad_COBS);
ra_ID_y = magFac*(secondarradiusYY + pad_COBS);
cx_ID = magFac*secondaryCenterX;
cy_ID = magFac*secondaryCenterY;
cxy = rotMat*[cx_ID; cy_ID];
cx_ID = cxy(1) + xShear;% + cshift; --> cshift not needed because falco_gen_ellipse takes centering into account
cy_ID = cxy(2) + yShear;% + cshift; --> cshift not needed because falco_gen_ellipse takes centering into account
% bm = prop_elliptical_obscuration(bm, ra_ID_x, ra_ID_y,'XC',cx_ID,'YC',cy_ID);

inputsSec.Narray = Narray;
inputsSec.Nbeam = Nbeam;
inputsSec.radiusX = ra_ID_x;
inputsSec.radiusY = ra_ID_y;
inputsSec.clockingDegrees = clock_deg;
inputsSec.centering = centering;
inputsSec.xShear = cx_ID;
inputsSec.yShear = cy_ID;
% inputsSec.magFac = magFac; %  not needed because falco_gen_ellipse takes magnification into account
secondaryObscuration = 1 - falco_gen_ellipse(inputsSec);

%% Struts

for iStrut=1:6
    angDeg = angStrutVec(iStrut) + clock_deg; % degrees
    wStrut = magFac*(wStrutVec(iStrut)+2*pad_strut);
    lStrutIn = magFac*lStrut;
    xc = magFac*(xcStrutVec(iStrut)); 
    yc = magFac*(ycStrutVec(iStrut)); 
    cxy = rotMat*[xc; yc];
    xc = cxy(1)+xShear;
    yc = cxy(2)+yShear;
    bm = prop_rectangular_obscuration(bm, lStrutIn, wStrut, 'XC',xc+cshift, 'YC',yc+cshift, 'ROTATION',angDeg);%, norm)
end

%% Tabs where Struts Meet Secondary Mirror

nTabs = 3;
tabCube = ones(Narray, Narray, nTabs);

for iTab = 1:nTabs
    cobsTabsMask = zeros(Narray);

    XSnew = (XS + tabCenterVecX(iTab)) - xShear;
    YSnew = (YS + tabCenterVecY(iTab)) - yShear;
    THETAS = atan2(YSnew,XSnew);
    
    if(angTabStart(iTab)>angTabEnd(iTab))
        cobsTabsMask( THETAS>=angTabEnd(iTab)+clock_rad & THETAS<=angTabStart(iTab)+clock_rad ) = 1.0;
    else
        cobsTabsMask( THETAS<=angTabEnd(iTab)+clock_rad & THETAS>=angTabStart(iTab)+clock_rad ) = 1.0;
    end

    %--ELLIPSE:
%     bm2 = prop_begin(Dbeam, wl, Narray,'beam_diam_fraction',bdf);

    % Full ellipse to be multiplied by the mask to get just tabs
    cx_tab = magFac*tabCenterVecX(iTab);
    cy_tab = magFac*tabCenterVecY(iTab);
    cxy = rotMat*[cx_tab; cy_tab];
    cx_tab = cxy(1)+xShear;
    cy_tab = cxy(2)+yShear;
    tabRadiusX = magFac*(tabRadiusVecX(iTab) + pad_COBStabs);
    tabRadiusY = magFac*(tabRadiusVecY(iTab) + pad_COBStabs);

%     bm2 = prop_elliptical_obscuration(bm2, tabRadiusX, tabRadiusY,'XC',cx_tab+cshift,'YC',cy_tab+cshift);
%     tabEllipse = 1-ifftshift(abs(bm2.wf));

    clear inputs
    inputs.Narray = Narray;
    inputs.Nbeam = Nbeam;% = 1219.0/ODpixels %Nbeam;
    inputs.radiusX = tabRadiusX;
    inputs.radiusY = tabRadiusY;
    inputs.clockingDegrees = clock_deg;
    inputs.centering = centering;
    inputs.xShear = cx_tab;
    inputs.yShear = cy_tab;
    % inputs.magFac = magFac;
    tabEllipse = falco_gen_ellipse(inputs);

    tabSector = cobsTabsMask.*tabEllipse;

    tabCube(:,:,iTab) = 1-tabSector;   
%     figure(11); imagesc(cobsTabsMask); axis xy equal tight; colorbar; drawnow;
%     figure(12); imagesc(tabEllipse); axis xy equal tight; colorbar; drawnow;
%     figure(13); imagesc(tabCube(:,:,iTab)); axis xy equal tight; colorbar; drawnow;

end

%% Output

pupil = secondaryObscuration.*primaryAperture.*tabCube(:,:,1).*tabCube(:,:,2).*tabCube(:,:,3).*ifftshift(abs(bm.wf));

if(flagRot180)
   pupil = rot90(pupil,2); 
end

end % EOF