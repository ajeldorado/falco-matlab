% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function pupil = falco_gen_pupil_LUVOIR_A_0(inputs,varargin)
%
%--Function to generate the LUVOIR Design A (#0 for the initial design) 
% telescope pupil in Matlab using PROPER.
%
% Coordinates and dimensions of the struts, primary, secondary, and 
% hex segments are from Matthew Bolcar at NASA GSFC:
%
% The aperture is a six-ring hexagonal array, with the inner most ring removed.  
% There are 120 segments.  Segments are 1.15-m flat-to-flat and segment gaps are 6 mm.  
% The struts are each 125 mm wide. The three struts are aligned vertically.  
% The bottom two struts are separated by 2.69538 m, center-to-center.
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Modified on 2018-05-04 by A.J. Riggs to be a function. 
% Written on 2017-09-07 by A.J. Riggs to generate the LUVOIR pupil. 

function pupil = falco_gen_pupil_LUVOIR_A_0(inputs,varargin)

%-------------------
% %--FOR DEBUGGING ONLY
% clear all;
% addpath ~/Repos/FALCO/lib/PROPER/
% inputs.Nbeam = 2000;
% inputs.centering = 'pixel';
% %-------------------

%--USER INPUTS
Nbeam   = inputs.Nbeam; % number of points across FULL usable pupil
centering = inputs.centering;% 'pixel' or 'interpixel' centering of the array

%--Primary mirror properties
width_hex = 1.15;
nrings = 6;
hexrad = 2/sqrt(3)*width_hex/2;
width_gap = 6e-3; % gap size between segments
width_strut = 125e-3; % meters
hexsep = width_hex + width_gap; 
Dap = 13*width_hex + 12*width_gap;
dx = Dap/Nbeam;

if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam+1/2); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even(Nbeam); %--number of points across output array. Same size as width when interpixel centered.
end
Darray = Narray*dx;

%--Values for PROPER
wl_dummy   = 1e-6;               % wavelength (m)
bdf = 1; %--beam diameter factor in output array

%--Centering of the aperture on the array
switch centering % 0 for pixel-centered pupil, or -dx/2 for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
    otherwise
        error('falco_gen_pupil_LUVOIR_A_0.m: Error! Centering must be either pixel or interpixel.')
end

%-------- Generate the input pupil for LUVOIR with PROPER
bm = prop_begin(Dap, wl_dummy, Narray,'beam_diam_fraction',bdf);

% Subtract the inner ring from all the rings
[bm,ap] = prop_hex_wavefront(bm,nrings,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
[~,ap2] = prop_hex_wavefront(bm,1,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
bm.wf = fftshift(ap-ap2);

%--Add the struts
bm = prop_rectangular_obscuration(bm, width_strut, 8*width_hex, 'XC',cshift, 'YC',cshift + Dap/4);
bm = prop_rectangular_obscuration(bm, width_strut, 8*width_hex, 'XC',cshift - 2.69538/2, 'YC',cshift - Dap/4);
bm = prop_rectangular_obscuration(bm, width_strut, 8*width_hex, 'XC',cshift + 2.69538/2, 'YC',cshift - Dap/4);

pupil = fftshift(bm.wf);
% figure(4); imagesc(pupil); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;


end %---END OF FUNCTION


% %--DEBUGGING: Visually verify that mask is centered correctly
% mask = pupil;
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

