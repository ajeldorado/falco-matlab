% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate the LUVOIR Design A, Aperture 5, telescope pupil in
%   Matlab using PROPER
% Coordinates and dimensions of the primary, secondary, and hex segments
%   are from Matthew Bolcar (NASA GSFC).
% Coordinates and dimensions of the secondary mirror support struts were a
%   best-fit match by A.J. Riggs by matching PROPER-made rectangles to the 
%   pupil file from Matthew Bolcar (NASA GSFC).
%
%--Coordinates of hex segments to skip:
% 1 13 114 115 126 127
% 1 12 113 114 125 126
%

function mask = falco_gen_pupil_LUVOIR_A_5_Lyot_struts(inputs,varargin)
% 
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
            error('falco_gen_pupil_LUVOIR_A_5_Lyot_struts: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

%--Centering of array: 'pixel' or 'interpixel'
if(isfield(inputs,'centering'))
    centering = inputs.centering;
else %--Default to pixel centering
    centering = 'pixel';
end

Dbeam = inputs.Dbeam; %--diameter of the beam at the mask (meters)
Nbeam   = inputs.Nbeam;     % number of points across the incoming beam           
ID = inputs.ID; % inner diameter of mask (in pupil diameters)
OD = inputs.OD; % outer diameter of mask (in pupil diameters)
wStrut = inputs.wStrut; % width of the struts (in pupil diameters)
wStrut = wStrut*Dbeam; %--now in meters
dx = Dbeam/Nbeam;

%--USER INPUTS
Nap   = Nbeam;%250;%324;%1500;%250;                  % number of points across FULL usable pupil

width_hex0 = 1.242; %-- flat-to-flat (m)
hexgap0 = 6e-3; % (m)
Dap = (12*width_hex0 + 12*hexgap0);
dx_drawing = 1.242/158; % (m) %--In actual drawing, 158 pixels across the 1.242m, so dx_pixel =
%dx_drawing. strut in drawing is 19 pixels across, so dx_strut =
%19*dx_drawing = 0.1494m

%--Other INPUTS
clock_deg = 0;%inputs.clock_deg; % clocking angle of the pupil (in degrees)
magfacD = 1;%inputs.magfacD; %magnification factor of the pupil diameter
xshift = 0;%inputs.xshift; % translation in x of pupil (in diameters)
yshift = 0;%inputs.yshift; % translation in y of pupil (in diameters)
pad_strut = 0; %2*pad_strut_pct/100*diam; %--Convert to meters. Factor of 2x needed since strut is padded on both sides

Dmask = Dbeam; % width of the beam (so can have zero padding if LS is undersized) (meters)

if(strcmpi(centering,'pixel') || strcmpi(centering,'odd'))
    Narray = ceil_even(Nbeam/magfacD + 1); % minimum even number of points across to fully contain the actual aperture (if interpixel centered)
else
    Narray = ceil_even(Nbeam/magfacD); % minimum even number of points across to fully contain the actual aperture (if interpixel centered)
end
Darray = Narray*dx; %--width of the output array (meters)

switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
        if(flagRot180deg)
            cshift = -dx;
        end
end

%--For PROPER 
wl_dummy   = 1e-6;               % wavelength (m)
bdf = Nbeam/Narray; %--beam diameter factor in output array
dx_t = 0;
dy_t = 0;
magfac = magfacD;

width_hex = magfac*width_hex0; %-- flat-to-flat (m)
nrings = 6;
hexrad = 2/sqrt(3)*width_hex/2;
hexgap = magfac*hexgap0; % (m)
hexsep = width_hex + hexgap; % distance from center to center of neighboring segments

wStrut0 = 19*dx_drawing*magfac; %%%%125e-3; % meters
%-------- Generate the input pupil for LUVOIR
bm = prop_begin(Dbeam, wl_dummy, Narray,'beam_diam_fraction',bdf);

%--PRIMARY MIRROR (OUTER DIAMETER)
ra_OD = (Dbeam*OD/2)*magfacD;
cx_OD = 0 + cshift + xshift;
cy_OD = 0 + cshift + yshift;
bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);

%--SECONDARY MIRROR (INNER DIAMETER)
ra_ID = (Dbeam*ID/2)*magfacD;
cx_ID = 0 + cshift + xshift;
cy_ID = 0 + cshift + yshift;
bm = prop_circular_obscuration(bm, ra_ID,'cx',cx_ID,'cy',cy_ID);

bm2 = bm;

% %--Add the struts

if(wStrut>0)
    bm2 = prop_rectangular_obscuration(bm2, wStrut, 7*width_hex, 'XC',cshift-dx_t, 'YC',cshift-dy_t + magfac*Dap/4);

    len_1a = 2*width_hex - 12*dx_drawing;
    bm2 = prop_rectangular_obscuration(bm2, wStrut, len_1a, 'XC',cshift-dx_t + (hexrad-0.5*wStrut0), 'YC',cshift-dy_t - len_1a/2.);
    bm2 = prop_rectangular_obscuration(bm2, wStrut, len_1a, 'XC',cshift-dx_t - (hexrad-0.5*wStrut0), 'YC',cshift-dy_t - len_1a/2.);

    len_1b = 3.75*width_hex;
    bm2 = prop_rectangular_obscuration(bm2, wStrut, len_1b, 'XC',cshift-dx_t + 1.25*hexrad*2, 'YC',cshift-dy_t - 3.5*width_hex,'ROT',30);
    bm2 = prop_rectangular_obscuration(bm2, wStrut, len_1b, 'XC',cshift-dx_t - 1.25*hexrad*2, 'YC',cshift-dy_t - 3.5*width_hex,'ROT',-30);
end

mask = ifftshift(abs(bm2.wf));
if(flagRot180deg)
    mask = rot90(mask,2);
end

end %--END OF FUNCTION
