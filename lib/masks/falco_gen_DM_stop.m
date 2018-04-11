% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate a circular aperture to place centered on the beam at a deformable mirror.
%
% Created by A.J. Riggs on 2017-11-15.
%
% INPUTS: 
%  rEdgesLeft:  vector of leading (i.e., rising) radial coordinate values of each transmissive ring
%  rEdgesRight: vector of trailing (i.e., falling) radial coordinate values of each transmissive ring
%  dx:      spatial resolution for a pixel (any units, but must be same as Dstop)
%  Dstop:     diameter of the shaped pupil (or rather of the beam assumed by the shaped pupil)
%  centering:   centering of beam in array. Either 'pixel' or 'interpixel'
%  (optional, to be added later if desired) inputs.flagBinary: flag to specify if pupil should be binary
%
% OUTPUTS:
%  DMstop:     2-D square array of a circular stop at a DM. Cropped down to the smallest even-sized array with no extra zero padding. 

function DMstop = falco_gen_DM_stop(dx,Dmask,centering)

% %--DEBUGGING ONLY: HARD-CODED INPUTS
% clear; close all;
% dx = 1e-3;
% Dstop = 48e-3;
% centering = 'interpixel';
% addpath ~/Repos/FALCO/proper_v3.0.1_matlab_22aug17/

diam = Dmask;% diameter of the mask (meters)
NapAcross = Dmask/dx; % minimum even number of points across to fully contain the actual aperture (if interpixel centered)
if(strcmpi(centering,'pixel'))
    Narray = 2*ceil(1/2*(Dmask/dx+1/2)); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = 2*ceil(1/2*Dmask/dx); %--number of points across output array. Same size as width when interpixel centered.
end

Darray = Narray*dx; %--width of the output array (meters)
bdf = 1; %--beam diameter factor in output array
wl_dummy   = 1e-6;     % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -Darray/2/Narray; 
    case {'pixel'}
        cshift = 0;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%--INITIALIZE PROPER
bm = prop_begin(Darray, wl_dummy, Narray,'beam_diam_fraction',bdf);
% figure(1); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;

%--Outer diameter of aperture
ra_OD = (diam/2); 
cx_OD = 0 + cshift;
cy_OD = 0 + cshift;

bm = prop_circular_aperture(bm, ra_OD,'XC',cx_OD,'YC',cy_OD);%, cx, cy, norm);
% figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;
% figure(3); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

DMstop = ifftshift(abs(bm.wf));
% figure(10); imagesc(DMstop); axis xy equal tight; colorbar; drawnow;
% figure(11); imagesc(DMstop-rot90(DMstop,2)); axis xy equal tight; colorbar; drawnow;


end %--END OF FUNCTION
