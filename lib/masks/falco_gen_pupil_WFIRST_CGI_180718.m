% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% The research was carried out at the Jet Propulsion Laboratory, 
% California Institute of Technology, under a contract with the 
% National Aeronautics and Space Administration.
% -------------------------------------------------------------------------
%
%  Function to generate WFIRST pupil 180718 in Matlab using PROPER 
%    functions for rectangles and circles.
%
%  Parameters were found by fitting circles and rectangles to the 
%    universally released file 'CGI_Entrance_Pupil_180718_No-labels.bmp'.
%
% %--Order of Operations:
% 1) Pad
% 2) Magnify
% 3) Rotate
% 4) Translate
%
% NOTE: All pupil features have normalized length units of pupil diameters.
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Created on 2018-08-07 by A.J. Riggs, Jet Propulsion Laboratory, California Institute of Technology.

function pupil = falco_gen_pupil_WFIRST_CGI_180718(Nbeam,centering,varargin)


%% Define the best-fit values for circles and rectangles (DO NOT CHANGE)
%--Values were found by fitting to the full-resolution, 4096x4096 pupil
%file, 'CGI_Entrance_Pupil_180718_No-labels.bmp'.

OD = 1.000130208333333;
xcOD = 8.680555555555557e-06;
ycOD = 8.680555555555557e-06;
ID = 3.030133333333332e-01;
xcCOBS = -1.155555555555556e-04;
ycCOBS = -6.133333333333334e-04;
IDtabs = 3.144078947368421e-01;
xcCOBStabs = -1.973684210526340e-04;
ycCOBStabs = -6.250000000000000e-03;
wStrutVec = [...
     3.219259259259259e-02; ...
     3.219259259259259e-02; ...
     3.219259259259259e-02; ...
     3.219259259259258e-02; ...
     3.219259259259259e-02; ...
     3.219259259259259e-02; ...
     ];
lStrut = 5.500000000000000e-01;
angStrutVec = [...
     4.308638741879215e+01; ...
     1.828091850580443e+01; ...
    -7.736372240624411e+01; ...
     7.746228722667239e+01; ...
    -1.833049685311381e+01; ...
    -4.310697246349373e+01; ...
    ];
xcStrutVec = [...
     1.637164789600492e-01; ...
     3.311169704392094e-01; ...
     1.542050924925356e-01; ...
    -1.556442459316893e-01; ...
    -3.075636241385107e-01; ...
    -1.712399202747162e-01; ...
    ];
ycStrutVec = [...
     2.695837795868052e-01; ...
     7.744558909460633e-03; ...
    -2.885875977555251e-01; ...
    -2.874651682155463e-01; ...
    -7.319997758726773e-04; ...
     2.748434070552074e-01; ...
     ];
angTabStart = [...
     1.815774989921760e+00; ...
    -3.487710035839058e-01; ...
    -2.416523875732038e+00; ...
    ];
angTabEnd = [...
     1.344727938801013e+00; ...
    -7.527300509955320e-01; ...
    -2.822938064533701e+00; ...
    ];

%% Changes to the pupil

if(size(varargin,2)>1)
    error('falco_gen_pupil_WFIRST_CGI_180718.m: Too many inputs')
elseif(size(varargin,2)==1)
    changes = varargin{1}; %--Structure containing which values to change;
else
    changes.dummy = 1; %--Else initialize structure
end
 
%% Oversized strut features: overwrite defaults if values specified.

if(isfield(changes,'OD')); OD = changes.OD;  end
if(isfield(changes,'ID')); ID = changes.ID;  end
if(isfield(changes,'wStrut')); wStrutVec = changes.wStrut*ones(6,1);  end
if(isfield(changes,'wStrutVec')); wStrutVec = changes.wStrutVec;  end %--Unlikely to be called. Would overrule line above this one.

%% Padding values for obscurations
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
flagRot180 = changes.flagRot180;
% % if(flagRot180)
% %     clock_deg = clock_deg-180;
% % end

%--Padding values. (pad_all is added to all the rest)
pad_all = changes.pad_all;%0.2/100; %--Uniform padding on all features
pad_strut = changes.pad_strut + pad_all;
pad_COBS = changes.pad_COBS + pad_all;
pad_COBStabs = changes.pad_COBStabs + pad_all;
pad_OD = changes.pad_OD + pad_all; %--Radial padding at the edge
 
% %--Padding values:
% pad_all = 0.2/100; %--Uniform padding on all features
% pad_strut = 2.9e-3/Dtel + pad_all;
% pad_COBS = 6e-3/Dtel + pad_all;
% pad_COBStabs = 6e-3/Dtel + pad_all;
% pad_OD = 0.3/100 + pad_all; %--Radial padding at the edge

%--Rotation matrix used on center coordinates.
rotMat = [cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)];

%% Coordinates

if(strcmpi(centering,'pixel') || strcmpi(centering,'odd'))
    Narray = ceil_even(Nbeam + 1); %--number of points across output array. Requires two more pixels when pixel centered.
else
    Narray = ceil_even(Nbeam); %--No zero-padding needed if beam is centered between pixels
end

if(strcmpi(centering,'interpixel'))
    xs = (-(Nbeam-1)/2:(Nbeam-1)/2)/Nbeam;
else
    xs = (-(Narray/2):Narray/2-1)/Nbeam; 
end

[XS,YS] = meshgrid(xs);
% RS = sqrt(XS.^2 + YS.^2);
% THETAS = atan2(YS,XS);

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

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%% INITIALIZE PROPER
bm = prop_begin(Dbeam, wl, Narray,'beam_diam_fraction',bdf);
% figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;

%% PRIMARY MIRROR (OUTER DIAMETER)
ra_OD = magFac*(OD/2-pad_OD);
cx_OD = magFac*xcOD;
cy_OD = magFac*ycOD;
cxy = rotMat*[cx_OD; cy_OD];
cx_OD = cxy(1)-xShear;
cy_OD = cxy(2)-yShear;
bm = prop_circular_aperture(bm, ra_OD,'XC',cx_OD+cshift,'YC',cy_OD+cshift);
% figure(3); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

%% SECONDARY MIRROR (INNER DIAMETER)
ra_ID = magFac*(ID/2 + pad_COBS);
cx_ID = magFac*(xcCOBS);
cy_ID = magFac*(ycCOBS);
cxy = rotMat*[cx_ID; cy_ID];
cx_ID = cxy(1)-xShear;
cy_ID = cxy(2)-yShear;
bm = prop_circular_obscuration(bm, ra_ID,'XC',cx_ID+cshift,'YC',cy_ID+cshift);
% figure(4); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

%% Struts

for istrut=1:6
    angDeg = angStrutVec(istrut) + clock_deg; % degrees
    wStrut = magFac*(wStrutVec(istrut)+pad_strut);
    lStrutIn = magFac*lStrut;
    xc = magFac*(xcStrutVec(istrut)); 
    yc = magFac*(ycStrutVec(istrut)); 
    cxy = rotMat*[xc; yc];
    xc = cxy(1)-xShear;
    yc = cxy(2)-yShear;
    bm = prop_rectangular_obscuration(bm, lStrutIn, wStrut, 'XC',xc+cshift, 'YC',yc+cshift, 'ROTATION',angDeg);%, norm)
    % figure(5); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;
end

%% TABS ON SECONDARY MIRROR
%--Compute as new shape, and then multiply the obscuration with the rest of
%the pupil.



%--SOFTWARE MASK:
XSnew = (1/1*XS+xcCOBStabs)+xShear;
YSnew = (1/1*YS+ycCOBStabs)+yShear;
% RSnew = (XSnew).^2 + (YSnew).^2;

overSizeFac = 1.3;
cobsTabsMask = zeros(Narray);
THETAS = atan2(YSnew,XSnew);
clock_rad = deg2rad(clock_deg);

if(angTabStart(1)>angTabEnd(1))
    cobsTabsMask( (XSnew).^2 + (YSnew).^2 <= (overSizeFac*magFac*IDtabs/2)^2 & ...
        ( ( THETAS>=angTabEnd(1)+clock_rad & THETAS<=angTabStart(1)+clock_rad )...
        | ( THETAS>=angTabEnd(2)+clock_rad & THETAS<=angTabStart(2)+clock_rad )...
        | ( THETAS>=angTabEnd(3)+clock_rad & THETAS<=angTabStart(3)+clock_rad ) )...
        ) = 1;
else
    cobsTabsMask( (XSnew).^2 + (YSnew).^2 <= (overSizeFac*magFac*IDtabs/2)^2 & ...
        ( ( THETAS<=angTabEnd(1)+clock_rad & THETAS>=angTabStart(1)+clock_rad )...
        | ( THETAS<=angTabEnd(2)+clock_rad & THETAS>=angTabStart(2)+clock_rad )...
        | ( THETAS<=angTabEnd(3)+clock_rad & THETAS>=angTabStart(3)+clock_rad ) )...
        ) = 1;
end


%--CIRCLE:
%--Initialize PROPER
bm2 = prop_begin(Dbeam, wl, Narray,'beam_diam_fraction',bdf);

%--Full circle of COBS tabs--to be multiplied by the mask to get just tabs
ra_tabs = magFac*(IDtabs/2 + pad_COBStabs);
cx_tabs = magFac*(xcCOBStabs);
cy_tabs = magFac*(ycCOBStabs);
cxy = rotMat*[cx_tabs; cy_tabs];
cx_tabs = cxy(1)-xShear;
cy_tabs = cxy(2)-yShear;
bm2 = prop_circular_obscuration(bm2, ra_tabs,'XC',cx_tabs+cshift,'YC',cy_tabs+cshift);

temp = 1-ifftshift(abs(bm2.wf));
temp = cobsTabsMask.*temp;

cobsTabs = 1-temp;
% figure(111); imagesc(cobsTabs); axis xy equal tight; colorbar; drawnow;
% figure(112); imagesc(RSnew); axis xy equal tight; colorbar; drawnow;
% figure(113); imagesc(THETAS); axis xy equal tight; colorbar; drawnow;


%% Output

pupil = cobsTabs.*ifftshift(abs(bm.wf));

if(flagRot180)
   pupil = rot90(pupil,2); 
end

% figure(11); imagesc(pupil); axis xy equal tight; colorbar; drawnow;


% %--Write out pupil to file.
% fn_out = 
% dlmwrite(fn_out,pupil,'precision','%d','delimiter',' ');
% pause(0.1)

end % EOF









