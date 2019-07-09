% Copyright 2018,2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% Function to generate the Keck telescope pupil using PROPER functions.
%
% Pupil and Lyot stop using the normalized distances (in telescope pupil diameters)
%  -Refer to Huby et al. 2017 (1701.06397.pdf) for list of pupil parameters 
%   in normalized pupil units.
%
% INPUT:
% - inputs = a structure with required field 'Nbeam' to set the number of
% points across the pupil diameter (flat-to-flat). Optional field
% 'centering' to set the pixel centering.
%
% OPTIONAL INPUT:
% - 'rotation',ang = pair of values making the pupil rotated by ang degrees. 
%
% OUTPUT:
% - pupil = the 2-D pupil amplitude map
%
% Updated on 2019-06-21 by A.J. Riggs.
% Created on 2017-05-09 by A.J. Riggs.

function pupil = falco_gen_pupil_Keck(inputs,varargin)

%--Set default values of input parameters
angRotDeg = 0;
%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'rotation'}
            icav = icav + 1;
            angRotDeg = angRotDeg+varargin{icav}; % For even arrays, beam center is in between pixels.
        otherwise
            error('falco_gen_pupil_Keck: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

%--Centering of array: 'pixel' or 'interpixel'
if(isfield(inputs,'centering'))
    centering = inputs.centering;
else %--Default to pixel centering
    centering = 'pixel';
end

% %%
% %-----------------------------------
% %--FOR DEBUGGING ONLY
% clear all; close all; clc;
% addpath ~/Documents/MATLAB/PROPER/
% angRotDeg = 4.5; %0;
% inputs.Nbeam = 500;
% % inputs.magfacD = 1;
% centering = 'pixel';
% %-----------------------------------

Dtel = 10.911920; % diameter of the primary used for sizing based on actual info (m)
diam = 1;  %--Diameter used for PROPER               
wl   = 3.776e-6;   % dummy value for wavelength (m)
np   = 2*inputs.Nbeam;                  % number of points

nrings=3;
hexrad=1.8/2/Dtel;
hexsep=hexrad*cos(pi/6)*2;
hexrad=hexrad-0.003/Dtel;
% nhex=nrings*(nrings+1)*3+1;

dx = diam/np; 
switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel','even'}
        cshift = -dx/2; 
    case {'pixel','odd'}
        cshift = 0;
%         if(flagRot180deg)
%             cshift = -dx;
%         end
end

strut_width_Pent = 0.0254/Dtel; % meters
ID_tel = 0.242854*diam; % meters

%-------- Generate the input pupil for Keck
bm = prop_begin(diam, wl, np);
% nh   = nrings * (nrings+1) * 3 + 1; % number of hexagonal segments

%--Make the hexes.
[bm] = prop_hex_wavefront(bm,nrings,hexrad,hexsep,'ROTATION',30+angRotDeg,'DARKCENTER','XCENTER',cshift,'YCENTER',cshift); %'NO_APPLY',

%--Make the secondary mirror obscuration
bm = prop_circular_obscuration(bm, ID_tel/2,'XC',cshift,'YC',cshift);

%--Make the secondary mirror support strut obscurations
bm = prop_rectangular_obscuration(bm, strut_width_Pent, diam, 'XC',cshift, 'YC',cshift, 'ROTATION',0+angRotDeg);
bm = prop_rectangular_obscuration(bm, strut_width_Pent, diam, 'XC',cshift, 'YC',cshift, 'ROTATION',60+angRotDeg);
bm = prop_rectangular_obscuration(bm, strut_width_Pent, diam, 'XC',cshift, 'YC',cshift, 'ROTATION',120+angRotDeg);

pupil = fftshift(bm.wf);
pupil = padOrCropEven(pupil,ceil_even(np/2*1.01)); %--Give a little extra padding in case the pupil is clocked by ~5 degrees so it doesn't exceed the array size.

% %--For DEBUGGING
% figure(2); imagesc(pupil); axis xy equal tight; title('Keck Pupil','Fontsize',20);
% figure(3); imagesc(pupil-rot90(pupil,2)); axis xy equal tight; colorbar; title('Input Pupil','Fontsize',20);
% figure(4); imagesc(padOrCropEven(pupil,np/2)); axis xy equal tight; title('Input Pupil','Fontsize',20);


end %--END OF FUNCTION