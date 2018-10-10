% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the LUVOIR Design A, Aperture 5, telescope pupil in
%   Matlab using PROPER
% Coordinates and dimensions of the primary, secondary, and hex segments
%   are from Matthew Bolcar (NASA GSFC).
% Coordinates and dimenstions of the secondary mirror support struts were a
%   best-fit match by A.J. Riggs by matching PROPER-made rectangles to the 
%   pupil file from Matthew Bolcar (NASA GSFC).
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Modified on 2018-02-25 by A.J. Riggs to be for LUVOIR A aperture 5. 
% Written on  2017-09-07 by A.J. Riggs to generate the first proposed LUVOIR pupil. 
%   Values for the geometry were provided by Matthew Bolcar at NASA GSFC.
%
%--Coordinates of hex segments to skip:
% 1 13 114 115 126 127
% 1 12 113 114 125 126


function mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs,varargin)

%--Set default values of input parameters
flagRot180deg = false;
%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'rot180'}
            flagRot180deg = true; % For even arrays, beam center is in between pixels.
        otherwise
            error('falco_gen_pupil_LUVOIR_A_5_mag_trans: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

% %-----------------------------------
% %--FOR DEBUGGING ONLY
% addpath ~/Repos/FALCO/lib/PROPER/
% inputs.Nbeam = 1000;
% inputs.magfacD = 1;
% flagRot180deg = true;
% inputs.centering = 'pixel';
% %-----------------------------------

%--Centering of array: 'pixel' or 'interpixel'
if(isfield(inputs,'centering'))
    centering = inputs.centering;
else %--Default to pixel centering
    centering = 'pixel';
end

%--Magnification factor of the pupil diameter
if(isfield(inputs,'magfacD')) 
    magfacD = inputs.magfacD;
else
    magfacD = 1;
end


%--Gap between primary mirror segments [meters]
if(isfield(inputs,'gap_width_m')) 
    hexgap0 = inputs.gap_width_m;
else
    hexgap0 = 6e-3; %--Default of 6.0 millimeters
end

Nbeam   = inputs.Nbeam;     % number of points across the incoming beam  




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% %




%--USER INPUTS
Nap   = Nbeam;%250;%324;%1500;%250;                  % number of points across FULL usable pupil
% magfac_vec = 1 + (-1:0.1:1)/100;
% magfac_vec = 1 + (-2:0.1:2)/100;


% Nmag = 1; %length(magfac_vec);


width_hex0 = 1.242; %-- flat-to-flat (m)
%hexgap0 = 6e-3; % (m)
Dap = (12*width_hex0 + 12*hexgap0);
dx = Dap/Nap;
dx_drawing = 1.242/158; % (m) %--In actual drawing, 158 pixels across the 1.242m, so dx_pixel =
%dx_drawing. strut in drawing is 19 pixels across, so dx_strut =
%19*dx_drawing = 0.1494m


if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam/magfacD+1/2); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even(Nbeam/magfacD); %--number of points across output array. Same size as width when interpixel centered.
end
Darray = Narray*dx;


%--For PROPER 
wl_dummy   = 1e-6;               % wavelength (m)
bdf = Nbeam/Narray; %--beam diameter factor in output array
      
dx_t = 0;%xyvec(ti,1);
dy_t = 0;%xyvec(ti,2);

width_hex = magfacD*width_hex0; %-- flat-to-flat (m)
nrings = 6;
hexrad = 2/sqrt(3)*width_hex/2;
hexgap = magfacD*hexgap0; % (m)
hexsep = width_hex + hexgap; % distance from center to center of neighboring segments


switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel','even'}
        cshift = -dx/2; 
    case {'pixel','odd'}
        cshift = 0;
        if(flagRot180deg)
            cshift = -dx; 
        end
end
strut_width0 = 19*dx_drawing*magfacD; %%%%125e-3; % meters
if(isfield(inputs,'strut_width'))
    strut_width = inputs.strut_width; % width of the struts (in pupil diameters)
    strut_width = strut_width*Darray; %--now in meters
else
    strut_width = strut_width0;
end




%-------- Generate the input pupil for LUVOIR
bm = prop_begin(Dap, wl_dummy, Narray,'beam_diam_fraction',bdf);

% Subtract the inner ring from all the rings
[ap] = falco_hex_aperture_LUVOIR_A_5(bm,nrings,hexrad,hexsep,'XC',cshift-dx_t,'YC',cshift-dy_t,'DARKCENTER'); %--Official Matlab PROPER from August 2017
% [~,ap2] = prop_hex_wavefront(bm,0,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
% bm.wf = fftshift(ap-ap2);

% %--Add the struts
bm = prop_rectangular_obscuration(bm, strut_width, 7*width_hex, 'XC',cshift-dx_t, 'YC',cshift-dy_t + magfacD*Dap/4);


len_1a = 2*width_hex - 12*dx_drawing;
bm = prop_rectangular_obscuration(bm, strut_width, len_1a, 'XC',cshift-dx_t + (hexrad-0.5*strut_width0), 'YC',cshift-dy_t - len_1a/2.);
bm = prop_rectangular_obscuration(bm, strut_width, len_1a, 'XC',cshift-dx_t - (hexrad-0.5*strut_width0), 'YC',cshift-dy_t - len_1a/2.);

len_1b = 3.75*width_hex;
bm = prop_rectangular_obscuration(bm, strut_width, len_1b, 'XC',cshift-dx_t + 1.25*hexrad*2, 'YC',cshift-dy_t - 3.5*width_hex,'ROT',30);
bm = prop_rectangular_obscuration(bm, strut_width, len_1b, 'XC',cshift-dx_t - 1.25*hexrad*2, 'YC',cshift-dy_t - 3.5*width_hex,'ROT',-30);


mask = ifftshift(abs(bm.wf)).*ap;

mask(mask>1) = 1; %--Get rid of overlapping segment edges at low resolution if the gap size is zero.
% figure(51); imagesc(mask); axis xy equal tight; colorbar;

if(flagRot180deg)
    mask = rot90(mask,2);
end


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


