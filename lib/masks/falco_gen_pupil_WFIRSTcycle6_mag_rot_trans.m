% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the Cycle 6 WFIRST CGI pupil (on-axis, though) in
% Matlab using PROPER.
% Coordinates and dimensions of the struts, primary, and secondary are from
% a Visio file from Kent Wallace of the March 2, 2016 Cycle 6 WFIRST pupil.
% The pupil used is on-axis.
%
% Modified by A.J. Riggs on 2017-02-09 to be a function with specifyable
%   amounts of magnification, clocking, and translation. 
%   Now keeping the pupil edges non-binary.
% This script was first written by A.J. Riggs on 2017-01-11.
%
%
%--Inputs:
% inputs.Dbeam: diameter of the incoming beam
% inputs.Npup: Number of points across the illuminated pupil diameter
% inputs.Narray: Number of points across the output array
% inputs.magfacD: magnification factor of the pupil diameter
% inputs.clock_deg: clocking angle of the pupil (in degrees)
% inputs.xshift: x-translation in (original) pupil diameters
% inputs.yshift: y-translation in (original) pupil diameters
%
%--To allow magnification, clocking, and translation:
% Center of each feature needs to be magnified and clocked first.
%   Then translate the centers.
% Circles can be magnified and translated in either order.
% Strut width also needs to scale with magnification


function mask = falco_gen_pupil_WFIRSTcycle6_mag_rot_trans(inputs,varargin)

% % %--DEBUGGING ONLY: HARD-CODED INPUTS
% clear
% inputs.Dbeam = 48e-3; % meters;
% inputs.Nbeam = 100;
% inputs.ID = 0.35;
% inputs.OD = 0.99;
% inputs.strut_width = 0.04;
% inputs.magfacD = 1;
% inputs.centering = 'pixel';
% inputs.xshift = 0;
% inputs.yshift = 0;
% inputs.clock_deg = 0;
% addpath(genpath('~/Repos/FALCO'))
% flagRot180deg = true;



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
            error('falco_gen_pupil_WFIRSTcycle6_mag_rot_trans: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

Dbeam = inputs.Dbeam; %--diameter of the beam at the mask (meters)
Nbeam   = inputs.Nbeam;     % number of points across the incoming beam  
% Narray = inputs.Narray;   % number of points across output array
clock_deg = inputs.clock_deg; % clocking angle of the pupil (in degrees)
magfacD = inputs.magfacD; %magnification factor of the pupil diameter
xshift = inputs.xshift; % translation in x of pupil (in diameters)
yshift = inputs.yshift; % translation in y of pupil (in diameters)

% Centering of array: 'pixel' or 'interpixel'
if(isfield(inputs,'centering'))
    centering = inputs.centering;
else %--Default to pixel centering
    centering = 'pixel';
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

pad_strut = 0;%2*pad_strut_pct/100*diam; %--Convert to meters. Factor of 2x needed
pad_idod = 0;%pad_idod_pct/100*diam; %--Convert to meters

dx = Dbeam/Nbeam;
Dmask = Dbeam; % % width of the mask (meters)
diam = Dmask;% width of the mask (meters)
% NapAcross = Dmask/dx; % minimum even number of points across to fully contain the actual aperture (if interpixel centered)
if(strcmpi(centering,'pixel'))
    Narray = 2*ceil(1/2*(Dmask/dx+1/2)); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = 2*ceil(1/2*Dmask/dx); %--number of points across output array. Same size as width when interpixel centered.
end

Darray = Narray*dx; %--width of the output array (meters)
bdf = 1; %--beam diameter factor in output array
wl_dummy   = 1e-6;     % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 shift for pixel-centered pupil, or -Darray/2/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -Darray/2/Narray; 
    case {'pixel'}
        cshift = 0;
        if(flagRot180deg)
            cshift = -Darray/Narray;
        end
end

%--DATA FROM THE VISIO FILE
D0 = 8; % inches, pupil diameter in Visio file
x0 = -26; % inches, pupil center in x in Visio file
y0 = 20.25; % inches, pupil center in y in Visio file
Dconv = (Dbeam/D0); % conversion factor from inches and Visio units to meters 

if(isfield(inputs,'strut_width'))
    strut_width = inputs.strut_width; %--strut width in pupil diameters
    strut_width = strut_width*Dbeam; %--now in meters
else
    strut_width = (0.209/D0)*Dbeam;
end


%--INITIALIZE PROPER
bm = prop_begin(Darray, wl_dummy, Narray,'beam_diam_fraction',bdf);
% figure(1); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;


%--PRIMARY MIRROR (OUTER DIAMETER)
ra_OD = (diam/2 - pad_idod)*magfacD;
cx_OD = 0 + cshift + xshift;
cy_OD = 0 + cshift + yshift;
% norm = false;
bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);%, cx, cy, norm);
% figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;
% figure(3); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


%--SECONDARY MIRROR (INNER DIAMETER)
ra_ID = ((2.46/2)*Dconv + pad_idod)*magfacD;
cx_ID = 0 + cshift + xshift;
cy_ID = 0 + cshift + yshift;
bm = prop_circular_obscuration(bm, ra_ID,'cx',cx_ID,'cy',cy_ID);%, cx, cy, norm)
% figure(4); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


% if(rot_deg0==0)
%     rot_coefs = 0;
% else
%     rot_coefs = [0,-1,1];
% end
% for ri = 1:length(rot_coefs)  %length(pad_rot_rad_vec);
%     rot_deg = rot_deg0*rot_coefs(ri);

    %--STRUT 1
    rot_s1 = 77.56 + clock_deg; % degrees
    sx_s1 = magfacD*((3.6/D0)*diam + pad_strut);
    sy_s1 = magfacD*(strut_width + pad_strut);
    cx_s1 = (-24.8566 - x0)*Dconv;
    cy_s1 = (22.2242 - y0)*Dconv;
    cxy1 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s1; cy_s1] - cshift) + cshift;
    cx_s1 = cxy1(1) + xshift + cshift;
    cy_s1 = cxy1(2) + yshift + cshift;
    bm = prop_rectangular_obscuration(bm, sx_s1, sy_s1, 'cx',cx_s1, 'cy',cy_s1,'rot', rot_s1);%, norm)
    % figure(5); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


    %--STRUT 2
    rot_s2 = -17.56 + clock_deg; % degrees
    sx_s2 = magfacD*((3.6/D0)*diam + pad_strut);
    sy_s2 = magfacD*(strut_width + pad_strut);
    cx_s2 = (-23.7187 - x0)*Dconv;
    cy_s2 = (20.2742 - y0)*Dconv;
    cxy2 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s2; cy_s2] - cshift) + cshift;
    cx_s2 = cxy2(1) + xshift + cshift;
    cy_s2 = cxy2(2) + yshift + cshift;
    bm = prop_rectangular_obscuration(bm, sx_s2, sy_s2, 'cx',cx_s2, 'cy',cy_s2, 'rot',rot_s2);%, norm)
    % figure(6); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


    %--STRUT 3
    rot_s3 = -42.44 + clock_deg; % degrees
    sx_s3 = magfacD*((3.6/D0)*diam + pad_strut);
    sy_s3 = magfacD*(strut_width + pad_strut);
    cx_s3 = (-24.8566 - x0)*Dconv;
    cy_s3 = (18.2758 - y0)*Dconv;
    cxy3 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s3; cy_s3] - cshift) + cshift;
    cx_s3 = cxy3(1) + xshift + cshift;
    cy_s3 = cxy3(2) + yshift + cshift;
    bm = prop_rectangular_obscuration(bm, sx_s3, sy_s3, 'cx',cx_s3, 'cy',cy_s3, 'rot',rot_s3);%, norm)
    % figure(7); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


    %--STRUT 4
    rot_s4 = 42.44 + clock_deg; % degrees
    sx_s4 = magfacD*((3.6/D0)*diam + pad_strut);
    sy_s4 = magfacD*(strut_width + pad_strut);
    cx_s4 = (-27.1434 - x0)*Dconv;
    cy_s4 = (18.2758 - y0)*Dconv;
    cxy4 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s4; cy_s4] - cshift) + cshift;
    cx_s4 = cxy4(1) + xshift + cshift;
    cy_s4 = cxy4(2) + yshift + cshift;
    bm = prop_rectangular_obscuration(bm, sx_s4, sy_s4, 'cx',cx_s4, 'cy',cy_s4, 'rot',rot_s4);%, norm)
    % figure(8); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

    %--STRUT 5
    rot_s5 = 17.56 + clock_deg; % degrees
    sx_s5 = magfacD*((3.6/D0)*diam + pad_strut);
    sy_s5 = magfacD*(strut_width + pad_strut);
    cx_s5 = (-28.2813 - x0)*Dconv;
    cy_s5 = (20.2742 - y0)*Dconv;
    cxy5 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s5; cy_s5] - cshift) + cshift;
    cx_s5 = cxy5(1) + xshift + cshift;
    cy_s5 = cxy5(2) + yshift + cshift;
    bm = prop_rectangular_obscuration(bm, sx_s5, sy_s5, 'cx',cx_s5, 'cy',cy_s5, 'rot',rot_s5);%, norm)
    % figure(9); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

    %--STRUT 6
    rot_s6 = 102.44 + clock_deg; % degrees
    sx_s6 = magfacD*((3.6/D0)*diam + pad_strut);
    sy_s6 = magfacD*(strut_width + pad_strut);
    cx_s6 = (-27.1434 - x0)*Dconv;
    cy_s6 = (22.2242 - y0)*Dconv;
    cxy6 = magfacD*[cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)]*([cx_s6; cy_s6] - cshift) + cshift;
    cx_s6 = cxy6(1) + xshift + cshift;
    cy_s6 = cxy6(2) + yshift + cshift;
    bm = prop_rectangular_obscuration(bm, sx_s6, sy_s6, 'cx',cx_s6, 'cy',cy_s6, 'rot',rot_s6);%, norm)
%     figure(10); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


mask = ifftshift(abs(bm.wf));
if(flagRot180deg)
    mask = rot90(mask,2);
end
% figure(19); imagesc(mask); axis xy equal tight; colorbar;


end %--END OF FUNCTION













