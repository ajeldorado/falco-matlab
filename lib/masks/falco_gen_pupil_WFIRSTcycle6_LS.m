% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the Cycle 6 WFIRST CGI pupil (on-axis, though) in
% Matlab using PROPER. This version allows padding of the ID, OD, and
% struts to use the mask as an under-sized Lyot stop.
% Coordinates and dimensions of the struts, primary, and secondary are from
% a Visio file from Kent Wallace of the March 2, 2016 Cycle 6 WFIRST pupil.
% The pupil used is on-axis.
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Modified by A.J. Riggs on 2017-02-09 to be a function with specifyable
%   amounts of magnification, clocking, and translation. 
%   Now keeping the pupil edges non-binary.
% This script was first written by A.J. Riggs on 2017-01-11.
%
%
%--Inputs:
% inputs.Dbeam: diameter of the incoming beam
% inputs.Nbeam: number of points across the incoming beam 
% inputs.ID: inner diameter of mask (in pupil diameters)
% inputs.OD: outer diameter of mask (in pupil diameters)
% inputs.strut_width: width of the struts (in pupil diameters)
% inputs.centering: centering of the pupil on the array ('pixel' or 'interpixel')
%
%--To allow magnification, clocking, and translation:
% Center of each feature needs to be magnified and clocked first.
%   Then translate the centers.
% Circles can be magnified and translated in either order.
% Strut width also needs to scale with magnification

function mask = falco_gen_pupil_WFIRSTcycle6_LS(inputs,varargin)

% Set default values of input parameters
flagRot180deg = false;
%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'rot180'}
            flagRot180deg = true; % For even arrays, beam center is in between pixels.
        otherwise
            error('falco_ gen_pupil_WFIRSTcycle6_LS: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end


% %--DEBUGGING ONLY: HARD-CODED INPUTS
% clear all
% addpath ~/Repos/FALCO/lib/PROPER/
% inputs.Dbeam = 46.3e-3; % meters;
% inputs.Nbeam = 124;
% inputs.ID = 0.40;
% inputs.OD = 0.80;
% inputs.strut_width = 0.04;
% inputs.centering = 'pixel';
% flagRot180deg = false;
% %--------------------------------------


Dbeam = inputs.Dbeam; %--diameter of the beam at the mask (meters)
Nbeam   = inputs.Nbeam;     % number of points across the incoming beam           
ID = inputs.ID; % inner diameter of mask (in pupil diameters)
OD = inputs.OD; % outer diameter of mask (in pupil diameters)
strut_width = inputs.strut_width; % width of the struts (in pupil diameters)

strut_width = strut_width*Dbeam; %--now in meters
dx = Dbeam/Nbeam;

% Centering of array: 'pixel' or 'interpixel'
if(isfield(inputs,'centering'))
    centering = inputs.centering;
else %--Default to pixel centering
    centering = 'pixel';
end

%--Other INPUTS
clock_deg = 0;%inputs.clock_deg; % clocking angle of the pupil (in degrees)
magfacD = 1;%inputs.magfacD; %magnification factor of the pupil diameter
xshift = 0;%inputs.xshift; % translation in x of pupil (in diameters)
yshift = 0;%inputs.yshift; % translation in y of pupil (in diameters)
pad_strut = 0; %2*pad_strut_pct/100*diam; %--Convert to meters. Factor of 2x needed since strut is padded on both sides


Dmask = Dbeam; % % width of the beam (so can have zero padding if LS is undersized) (meters)
diam = Dmask;% width of the mask (meters)


if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam+1/2); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even(Nbeam); %--number of points across output array. Same size as width when interpixel centered.
end
Darray = Narray*dx; %--width of the output array (meters)

bdf = Nbeam/Narray; %--beam diameter factor in output array
wl_dummy   = 1e-6;     % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
        if(flagRot180deg)
            cshift = -dx;
        end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%--DATA FROM THE VISIO FILE
D0 = 8; % inches, pupil diameter in Visio file
x0 = -26; % inches, pupil center in x in Visio file
y0 = 20.25; % inches, pupil center in y in Visio file
Dconv = (diam/D0); % conversion factor from inches and Visio units to meters 

%--Nominal strut width in units of normalized diameters.
% strut_width = 0.209*(diam/D0);

%--INITIALIZE PROPER
bm = prop_begin(Dbeam, wl_dummy, Narray,'beam_diam_fraction',bdf);
% figure(1); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;

%--PRIMARY MIRROR (OUTER DIAMETER)
% ra_OD = (diam/2 - pad_idod)*magfacD;
ra_OD = (Dbeam*OD/2)*magfacD;
cx_OD = 0 + cshift + xshift;
cy_OD = 0 + cshift + yshift;
% norm = false;
bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);%, cx, cy, norm);
% figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;
% figure(3); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

%--SECONDARY MIRROR (INNER DIAMETER)
ra_ID = (Dbeam*ID/2)*magfacD; %((2.46/2)*Dconv + pad_idod)*magfacD;
cx_ID = 0 + cshift + xshift;
cy_ID = 0 + cshift + yshift;
bm = prop_circular_obscuration(bm, ra_ID,'cx',cx_ID,'cy',cy_ID);%, cx, cy, norm)
% figure(4); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;
% figure(104); imagesc(ifftshift(abs(bm.wf))-rot90(ifftshift(abs(bm.wf)),2)); axis xy equal tight; colorbar;



% if(rot_deg0==0)
%     rot_coefs = 0;
% else
%     rot_coefs = [0,-1,1];
% end
% for ri = 1:length(rot_coefs)  %length(pad_rot_rad_vec);
%     rot_deg = rot_deg0*rot_coefs(ri);


    %--STRUT 1
    rot_s1 = 77.56 + clock_deg; % degrees
    sx_s1 = magfacD*(3.6*(diam/D0) + pad_strut);
    sy_s1 = magfacD*(strut_width + pad_strut);
    cx_s1 = (-24.8566 - x0)*Dconv + cshift;
    cy_s1 = (22.2242 - y0)*Dconv + cshift;
    cxy1 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s1; cy_s1] - cshift) + cshift;
    cx_s1 = cxy1(1) + xshift;
    cy_s1 = cxy1(2) + yshift;
    bm = prop_rectangular_obscuration(bm, sx_s1, sy_s1, 'cx',cx_s1, 'cy',cy_s1, 'rot', rot_s1);%, norm)
    % figure(5); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


    %--STRUT 2
    rot_s2 = -17.56 + clock_deg; % degrees
    sx_s2 = magfacD*(3.6*(diam/D0) + pad_strut);
    sy_s2 = magfacD*(strut_width + pad_strut);
    cx_s2 = (-23.7187 - x0)*Dconv + cshift;
    cy_s2 = (20.2742 - y0)*Dconv + cshift;
    cxy2 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s2; cy_s2] - cshift) + cshift;
    cx_s2 = cxy2(1) + xshift;
    cy_s2 = cxy2(2) + yshift;
    bm = prop_rectangular_obscuration(bm, sx_s2, sy_s2, 'cx',cx_s2, 'cy',cy_s2, 'rot', rot_s2);%, norm)
    % figure(6); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


    %--STRUT 3
    rot_s3 = -42.44 + clock_deg; % degrees
    sx_s3 = magfacD*(3.6*(diam/D0) + pad_strut);
    sy_s3 = magfacD*(strut_width + pad_strut);
    cx_s3 = (-24.8566 - x0)*Dconv + cshift;
    cy_s3 = (18.2758 - y0)*Dconv + cshift;
    cxy3 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s3; cy_s3] - cshift) + cshift;
    cx_s3 = cxy3(1) + xshift;
    cy_s3 = cxy3(2) + yshift;
    bm = prop_rectangular_obscuration(bm, sx_s3, sy_s3, 'cx',cx_s3, 'cy',cy_s3, 'rot', rot_s3);%, norm)
    % figure(7); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


    %--STRUT 4
    rot_s4 = 42.44 + clock_deg; % degrees
    sx_s4 = magfacD*(3.6*(diam/D0) + pad_strut);
    sy_s4 = magfacD*(strut_width + pad_strut);
    cx_s4 = (-27.1434 - x0)*Dconv + cshift;
    cy_s4 = (18.2758 - y0)*Dconv + cshift;
    cxy4 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s4; cy_s4] - cshift) + cshift;
    cx_s4 = cxy4(1) + xshift;
    cy_s4 = cxy4(2) + yshift;
    bm = prop_rectangular_obscuration(bm, sx_s4, sy_s4, 'cx',cx_s4, 'cy',cy_s4, 'rot', rot_s4);%, norm)
    % figure(8); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

    %--STRUT 5
    rot_s5 = 17.56 + clock_deg; % degrees
    sx_s5 = magfacD*(3.6*(diam/D0) + pad_strut);
    sy_s5 = magfacD*(strut_width + pad_strut);
    cx_s5 = (-28.2813 - x0)*Dconv + cshift;
    cy_s5 = (20.2742 - y0)*Dconv + cshift;
    cxy5 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s5; cy_s5] - cshift) + cshift;
    cx_s5 = cxy5(1) + xshift;
    cy_s5 = cxy5(2) + yshift;
    bm = prop_rectangular_obscuration(bm, sx_s5, sy_s5, 'cx',cx_s5, 'cy',cy_s5, 'rot', rot_s5);%, norm)
    % figure(9); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

    %--STRUT 6
    rot_s6 = 102.44 + clock_deg; % degrees
    sx_s6 = magfacD*(3.6*(diam/D0) + pad_strut);
    sy_s6 = magfacD*(strut_width + pad_strut);
    cx_s6 = (-27.1434 - x0)*Dconv + cshift;
    cy_s6 = (22.2242 - y0)*Dconv + cshift;
    cxy6 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s6; cy_s6] - cshift) + cshift;
    cx_s6 = cxy6(1) + xshift;
    cy_s6 = cxy6(2) + yshift;
    bm = prop_rectangular_obscuration(bm, sx_s6, sy_s6, 'cx',cx_s6, 'cy',cy_s6, 'rot', rot_s6);%, norm)
%     figure(10); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

%     pause(1)


mask = ifftshift(abs(bm.wf));
if(flagRot180deg)
    mask = rot90(mask,2);
end

% figure(18); imagesc(mask); axis xy equal tight; colorbar;
% figure(19); imagesc(mask-rot90(mask,2)); axis xy equal tight; colorbar;


end %--END OF FUNCTION




% %--DEBUGGING: Visually verify that mask is centered correctly
% figure(11); imagesc(mask); axis xy equal tight; colorbar; drawnow;
% switch centering 
%     case {'pixel'}
%         maskTemp = mask(2:end,2:end);
%     otherwise
%         maskTemp = mask;
% end
% figure(12); imagesc(maskTemp-fliplr(maskTemp)); axis xy equal tight; colorbar; 
% title('Centering Check','Fontsize',20); set(gca,'Fontsize',20);
% drawnow;
% 
% sum(sum(mask))








