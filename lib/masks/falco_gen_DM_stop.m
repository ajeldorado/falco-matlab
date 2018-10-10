% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate a circular aperture to place centered on the beam at a deformable mirror.
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Created on 2017-11-15 by A.J. Riggs (JPL).
%
% INPUTS: 
%  dx:      spatial resolution for a pixel [any units as long as they match that of Dmask]
%  Dmask:     diameter of the aperture mask [any units as long as they match that of dx]
%  centering:   centering of beam in array. Either 'pixel' or 'interpixel'
%
% OUTPUTS:
%  mask:     2-D square array of a circular stop at a DM. Cropped down to the smallest even-sized array with no extra zero padding. 

function mask = falco_gen_DM_stop(dx,Dmask,centering)

% %--DEBUGGING ONLY: HARD-CODED INPUTS
% clear; close all;
% dx = 1e-3;
% Dmask = 50e-3;%46.3e-3;
% centering = 'interpixel';
% addpath ~/Repos/FALCO/lib/PROPER/


Nbeam = Dmask/dx; %--Number of points across the mask.

%--Minimum number of points across the array to fully contain the mask
if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam+1/2); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even(Nbeam); %--number of points across output array. Same size as width when interpixel centered.
end

% Darray = Narray*dx; %--width of the output array (meters)
bdf = Nbeam/Narray; %--beam diameter factor in output array
wl_dummy   = 1e-6;     % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; % = -dx/2; 
    case {'pixel'}
        cshift = 0;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%--INITIALIZE PROPER. Note that:  bm.dx = diam / bdf / np;
bm = prop_begin(Dmask, wl_dummy, Narray,'beam_diam_fraction',bdf);
% figure(1); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;

%--Outer diameter of aperture
ra_OD = (Dmask/2); 
cx_OD = 0 + cshift;
cy_OD = 0 + cshift;

bm = prop_circular_aperture(bm, ra_OD,'XC',cx_OD,'YC',cy_OD);%, cx, cy, norm);
% figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;
% figure(3); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

mask = ifftshift(abs(bm.wf));

end %--END OF FUNCTION


% %--DEBUGGING: Visually verify that mask is centered correctly
% figure(11); imagesc(mask); axis xy equal tight; colorbar; drawnow;
% switch centering 
%     case {'pixel'}
%         maskTemp = mask(2:end,2:end);
%     otherwise
%         maskTemp = mask;
% end
% figure(12); imagesc(maskTemp-rot90(maskTemp,2)); axis xy equal tight; colorbar; 
% title('Centering Check','Fontsize',20); set(gca,'Fontsize',20);
% drawnow;
% 
% sum(sum(mask))