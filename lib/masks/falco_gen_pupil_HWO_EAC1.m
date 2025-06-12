% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate the HWO EAC1 telescope pupil in Matlab using PROPER.
%
% From another script:
% script_check_hexagonal_aperture_inscribed
% Inscribed diameter = 5972.13594623
% Circumscribed diameter = 7219.39682680

function pupil = falco_gen_pupil_HWO_EAC1(inputs)

% 4mm physical gap but 6mm optical gap, so take an extra 2mm off
% flat-to-flat mirror width.
flat_to_flat_physical_m = 1431.614*1e-3; % meters
segment_gap_physical_m = 4e-3; % meters
% flat_to_flat = (1431.614-2)*1e-3; % meters
% segment_gap_m = (4+2)*1e-3; % meters

if(isfield(inputs,'centering')); centering = inputs.centering; else; centering = 'pixel'; end % array centering
if(isfield(inputs,'magfacD')); mag = inputs.magfacD; else; mag = 1; end  %--Pupil Magnification
if(isfield(inputs,'segment_gap_m')); segment_gap_m = inputs.segment_gap_m; else; segment_gap_m = (4+2)*1e-3; end % meters
if(isfield(inputs,'clock_deg')); clockDeg = inputs.clock_deg; else; clockDeg = 0; end  % 23.4350 to make pupil edge touch array edge
if(isfield(inputs,'flagRot180deg')); flagRot180deg = inputs.flagRot180deg; else; flagRot180deg = false; end
if(isfield(inputs, 'upsampleFactor')); upsampleFactor = inputs.upsampleFactor; else; upsampleFactor = 100; end % upsampling factor
if(isfield(inputs,'xShear')); xShear = inputs.xShear; else; xShear = 0; end % [pupil diameters]
if(isfield(inputs,'yShear')); yShear = inputs.yShear; else; yShear = 0; end % [pupil diameters]
shearMax = max(abs([xShear, yShear])); % [pupil diameters]

% shrinnk the flat-to-flat by the same increase in segment gap so that the
% diameter stays the same.
flat_to_flat = flat_to_flat_physical_m - (segment_gap_m - segment_gap_physical_m);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--Fixed INPUTS
Nbeam0 = inputs.Nbeam;     % number of points across the incoming beam
scaleFac = 1; %0.96075;
Nbeam = Nbeam0;%*scaleFac;  % Change Nbeam to be flat-to-flat. Nbeam0 is for the circumscribing circle. makes the beam size match the hypergaussian approach
nrings = 2;
width_hex0 = flat_to_flat; %0.955; % flat-to-flat (m)
Dap = 7219.39682680*1e-3; % Circumscribed diameter in mm of HWO. % Only measured for 6mm gaps and (1431.614-2)mm flat to flat.
% Dap = ((2*nrings)*width_hex0 + (2*nrings-1)*segment_gap_m);
dx = Dap/Nbeam0;

xShearM = 1/scaleFac*Dap*xShear; % meters
yShearM = 1/scaleFac*Dap*yShear; % meters

if(strcmpi(centering,'pixel'))
    Narray = ceil_even((1+2*shearMax)*Nbeam*max([mag, 1])+1); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even((1+2*shearMax)*Nbeam*max([mag, 1])); %--number of points across output array. Same size as width when interpixel centered.
end
Narray = Narray+2;
% Darray = Narray*dx;

%--For PROPER 
wl_dummy = 1e-6; % wavelength (m)
bdf = Nbeam/Narray; %--beam diameter factor in output array

hexWidth = mag*width_hex0; %-- flat-to-flat (m)
hexRadius = 2/sqrt(3)*hexWidth/2;
hexGap = mag*segment_gap_m; % (m)
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
prop_set_antialiasing(ceil_odd(upsampleFactor));
[ap] = falco_hex_aperture_HWO_EAC1(bm, nrings, hexRadius, hexSep, 'XC', cShift+xShearM, 'YC', cShift+yShearM, 'ROTATION', clockDeg);

pupil = ifftshift(abs(bm.wf)).*ap;

pupil(pupil > 1) = 1; %--Get rid of overlapping segment edges at low resolution if the gap size is zero.

if(flagRot180deg)
    pupil = rot90(pupil, 2);
end

end %--END OF FUNCTION
