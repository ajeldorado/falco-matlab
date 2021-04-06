% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate the Keck telescope pupil using PROPER functions.
%
% INPUTS
% ------
% inputs : a structure of inputs
%  - inputs.Nbeam = number of points across the pupil diameter (flat-to-flat).
%  - inputs.centering : (optional) centering of the pupil in the array.
%  - inputs.clocking : (optional) CCW rotation angle of the pupil (degrees)
%
% OUTPUTS
% -------
% pupil : the 2-D pupil amplitude map
%
% REFERENCES
% ----------
% Huby et al. 2017 (1701.06397.pdf) provides the pupil specifications

function pupil = falco_gen_pupil_Keck(inputs)

%--Set default values of input parameters
clocking = 0; % [degrees CCW]
centering = 'pixel'; % 'interpixel' or 'pixel'
xShear = 0; % [pupil diameters]
yShear = 0; % [pupil diameters]

%--Centering of array: 'pixel' or 'interpixel'
if isfield(inputs, 'centering'); centering = inputs.centering; end % 'interpixel' or 'pixel'
if isfield(inputs, 'clocking'); clocking = inputs.clocking; end % [degrees CCW]
if isfield(inputs, 'xShear'); xShear = inputs.xShear; end % [pupil diameters]
if isfield(inputs, 'yShear'); yShear = inputs.yShear; end % [pupil diameters]

Dtel = 10.911920; % diameter of the primary used for sizing based on actual info (m)
diam = 1;  %--Diameter used for PROPER               
wl   = 3.776e-6;   % dummy value for wavelength (m)

mag = 1;
shearMax = max(abs([xShear, yShear])); % [pupil diameters]
if(strcmpi(centering,'pixel'))
    Narray = ceil_even(1.01*(1+2*shearMax)*inputs.Nbeam*max([mag, 1])+1); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even(1.01*(1+2*shearMax)*inputs.Nbeam*max([mag, 1])); %--number of points across output array. Same size as width when interpixel centered.
end

nrings=3;
hexrad=1.8/2/Dtel;
hexsep=hexrad*cos(pi/6)*2;
hexrad=hexrad-0.003/Dtel;
% nhex=nrings*(nrings+1)*3+1;

dx = diam/Narray; 
switch centering % 0 for pixel-centered pupil, or -diam/Narray for inter-pixel centering
    case {'interpixel','even'}
        cshift = -dx/2; 
    case {'pixel','odd'}
        cshift = 0;
end

strut_width_Pent = 0.0254/Dtel; % meters
ID_tel = 0.242854*diam; % meters

%-------- Generate the input pupil for Keck
bdf = inputs.Nbeam/Narray; %--beam diameter factor in output array

bm = prop_begin(diam, wl, Narray, bdf);
% nh   = nrings * (nrings+1) * 3 + 1; % number of hexagonal segments

%--Make the hexes.
bm = prop_hex_wavefront(bm,nrings,hexrad,hexsep,'ROTATION',30+clocking,'DARKCENTER','XCENTER',cshift+xShear,'YCENTER',cshift+yShear);

%--Make the secondary mirror obscuration
bm = prop_circular_obscuration(bm, ID_tel/2,'XC',cshift+xShear,'YC',cshift+yShear);

%--Make the secondary mirror support strut obscurations
bm = prop_rectangular_obscuration(bm, strut_width_Pent, diam, 'XC',cshift+xShear, 'YC',cshift+yShear, 'ROTATION',0+clocking);
bm = prop_rectangular_obscuration(bm, strut_width_Pent, diam, 'XC',cshift+xShear, 'YC',cshift+yShear, 'ROTATION',60+clocking);
bm = prop_rectangular_obscuration(bm, strut_width_Pent, diam, 'XC',cshift+xShear, 'YC',cshift+yShear, 'ROTATION',120+clocking);

pupil = fftshift(bm.wf);
% pupil = pad_crop(pupil, ceil_even(Narray/2*1.01)); %--Give a little extra padding so that the pupil never extends past the array bounds.

end %--END OF FUNCTION
