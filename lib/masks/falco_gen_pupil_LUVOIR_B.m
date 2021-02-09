% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate the LUVOIR Design B telescope pupil in Matlab using PROPER.
%
function pupil = falco_gen_pupil_LUVOIR_B(inputs,varargin)

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
if(isfield(inputs,'wGap_m')); hexGap0 = inputs.wGap_m; else; hexGap0 = 6e-3; end %--Default of 6.0 millimeters
if(isfield(inputs,'clock_deg')); clockDeg = inputs.clock_deg; else; clockDeg = 0; end
if(isfield(inputs,'xShear')); xShear = inputs.xShear; else; xShear = 0; end % [pupil diameters]
if(isfield(inputs,'yShear')); yShear = inputs.yShear; else; yShear = 0; end % [pupil diameters]
shearMax = max(abs([xShear, yShear])); % [pupil diameters]

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--USER INPUTS
Nbeam0 = inputs.Nbeam;     % number of points across the incoming beam
scaleFac = 0.96075;
Nbeam = Nbeam0*scaleFac;  % Change Nbeam to be flat-to-flat. Nbeam0 is for the circumscribing circle. makes the beam size match the hypergaussian approach
nrings = 4;
width_hex0 = 0.955; % flat-to-flat (m)
Dap = ((2*nrings)*width_hex0 + (2*nrings-1)*hexGap0);
dx = Dap/Nbeam0;

xShearM = 1/scaleFac*Dap*xShear; % meters
yShearM = 1/scaleFac*Dap*yShear; % meters

if(strcmpi(centering,'pixel'))
    Narray = ceil_even((1+2*shearMax)*1.02*Nbeam*max([mag, 1])+1); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even((1+2*shearMax)*1.02*Nbeam*max([mag, 1])); %--number of points across output array. Same size as width when interpixel centered.
end

% Darray = Narray*dx;

%--For PROPER 
wl_dummy = 1e-6; % wavelength (m)
bdf = Nbeam/Narray; %--beam diameter factor in output array

hexWidth = mag*width_hex0; %-- flat-to-flat (m)
hexRadius = 2/sqrt(3)*hexWidth/2;
hexGap = mag*hexGap0; % (m)
hexSep = hexWidth + hexGap; % distance from center to center of neighboring segments

switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel'}
        cShift = -dx/2; 
    case {'pixel'}
        cShift = 0;
        if(flagRot180deg)
            cShift = -dx;
        end
end

% Use PROPER to generate the hexagonal mirrors
bm = prop_begin(Dap, wl_dummy, Narray,'beam_diam_fraction',bdf);
[ap] = falco_hex_aperture_LUVOIR_B(bm, nrings, hexRadius, hexSep, 'XC', cShift+xShearM, 'YC', cShift+yShearM, 'ROTATION', clockDeg);

pupil = ifftshift(abs(bm.wf)).*ap;

pupil(pupil > 1) = 1; %--Get rid of overlapping segment edges at low resolution if the gap size is zero.

if(flagRot180deg)
    pupil = rot90(pupil, 2);
end

end %--END OF FUNCTION
