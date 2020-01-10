% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate the LUVOIR Design A (Final) telescope pupil from 
% 2018 in Matlab using PROPER.
%
% Coordinates and dimensions of the primary, secondary, and hex segments
%   are from Matthew Bolcar (NASA GSFC).
% Coordinates and dimensions of the secondary mirror support struts were a
%   best-fit match by A.J. Riggs by matching PROPER-made rectangles to the 
%   pupil file from Matthew Bolcar (NASA GSFC).
%
%--Coordinates of hex segments to skip:
% 1 13 114 115 126 127
% 1 12 113 114 125 126

function pupil = falco_gen_pupil_LUVOIR_A_final(inputs,varargin)

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

%--Centering of array: 'pixel' or 'interpixel'
if(isfield(inputs,'centering'))
    centering = inputs.centering;
else
    centering = 'pixel';
end

%--Pupil Magnification
if(isfield(inputs,'magfacD')) 
    mag = inputs.magfacD;
else
    mag = 1;
end

%--Gap between primary mirror segments [meters]
if(isfield(inputs,'wGap_m')) 
    hexGap0 = inputs.wGap_m;
else
    hexGap0 = 6e-3; %--Default of 6.0 millimeters
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--USER INPUTS
Nbeam   = inputs.Nbeam;     % number of points across the incoming beam  
Nap   = Nbeam; % number of points across FULL usable pupil
width_hex0 = 1.2225; %-- flat-to-flat (m)
Dap = (12*width_hex0 + 11*hexGap0);%(12*width_hex0 + 12*hexgap0);
dx = Dap/Nap;
dxDrawing = 1.2225/158; % (m) %--In actual drawing, 158 pixels across the 1.2225m, so dx_pixel =
%dx_drawing. strut in drawing is 19 pixels across, so dx_strut =
%19*dx_drawing = 0.1470m

if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam/mag+1); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even(Nbeam/mag); %--number of points across output array. Same size as width when interpixel centered.
end
Darray = Narray*dx;

%--For PROPER 
wl_dummy = 1e-6; % wavelength (m)
bdf = Nbeam/Narray; %--beam diameter factor in output array
      
dx_t = 0;
dy_t = 0;

hexWidth = mag*width_hex0; %-- flat-to-flat (m)
nrings = 6;
hexRadius = 2/sqrt(3)*hexWidth/2;
hexGap = mag*hexGap0; % (m)
hexSep = hexWidth + hexGap; % distance from center to center of neighboring segments

switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel','even'}
        cShift = -dx/2; 
    case {'pixel','odd'}
        cShift = 0;
        if(flagRot180deg)
            cShift = -dx;
        end
end
strutWidth = 0.15*mag; % meters

% Use PROPER to generate the hexagonal mirrors and rectangular struts.
bm = prop_begin(Dap, wl_dummy, Narray,'beam_diam_fraction',bdf);
[ap] = falco_hex_aperture_LUVOIR_A_5(bm,nrings,hexRadius,hexSep,'XC',cShift-dx_t,'YC',cShift-dy_t,'DARKCENTER');
bm = prop_rectangular_obscuration(bm, strutWidth, 7*hexWidth, 'XC',cShift-dx_t, 'YC',cShift-dy_t + mag*Dap/4);
len_1b = (sqrt(93)+0.5)*hexRadius;
bm = prop_rectangular_obscuration(bm, strutWidth, len_1b, 'XC',cShift-dx_t + 1.5*hexRadius, 'YC',cShift-dy_t - 11*sqrt(3)/4*hexRadius,'ROT',12.7);
bm = prop_rectangular_obscuration(bm, strutWidth, len_1b, 'XC',cShift-dx_t - 1.5*hexRadius, 'YC',cShift-dy_t - 11*sqrt(3)/4*hexRadius,'ROT',-12.7);

pupil = ifftshift(abs(bm.wf)).*ap;

pupil(pupil > 1) = 1; %--Get rid of overlapping segment edges at low resolution if the gap size is zero.

if(flagRot180deg)
    pupil = rot90(pupil, 2);
end

end %--END OF FUNCTION
