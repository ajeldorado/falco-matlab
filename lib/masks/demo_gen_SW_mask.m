% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate binary (0-1) software masks for the focal plane. 
% This can be used as a field stop, 
% or for making the scoring and correction regions in the focal plane.
%
% Created on 2018-03-07 by A.J. Riggs
%
%--INPUTS:
% inputs: structure with several fields:
%  -pixresFP: pixels per lambda_c/D
%  -rhoInner: radius of inner FPM amplitude spot (in lambda_c/D)
%  -rhoOuter: radius of outer opaque FPM ring (in lambda_c/D)
%  -angDeg: angular opening (degrees) on the left/right/both sides.
%  -whichSide: which sides to have open. 'left','right', 'top', 'bottom', or 'both'
%  -centering: centering of the coordinates. 'pixel' or 'interpixel'
%  -FOV: minimum desired field of view (in lambda_c/D). Omitting sets to rhoOuter
%  -shape: 'square' makes a square. Omitting makes a circle. 
%  -clockAngDeg: Dark hole rotation about the z-axis (deg)
%
%--OUTPUTS:
% maskSW: rectangular, even-sized, binary-valued software mask
% xis: vector of coordinates along the horizontal axis (in lambda_c/D)
% etas: : vector of coordinates along the vertical axis (in lambda_c/D)

clear all;

% pixresFP = 6;
% 
% centering = 'pixel';
% Nxi = 100;
% Neta = 100;
% xiOffset = 0;
% etaOffset = 0;
% 
% %--Focal Plane Coordinates
% dxi = 1/pixresFP;
% deta = dxi;
% if( strcmpi(centering,'interpixel')  )
%     xis  = (-(Nxi-1)/2: (Nxi-1)/2 )*dxi;
%     etas = (-(Neta-1)/2:(Neta-1)/2)*deta;
% else %--pixel centering
%     xis  = (-Nxi/2: (Nxi/2-1) )*dxi;
%     etas = (-Neta/2:(Neta/2-1))*deta;
% end
% [XIS, ETAS] = meshgrid(xis,etas);
% XIS = XIS - xiOffset;
% ETAS = ETAS - etaOffset;
% [THETAS, RHOS] = cart2pol(XIS, ETAS);
% 
% mask = RHOS < 6  & RHOS.*cos(THETAS) >= 1;
% 
% figure(11); imagesc(xis, etas, mask); axis xy equal tight; colormap parula; colorbar; drawnow;
% figure(12); imagesc(xis, etas, RHOS.*cos(THETAS)); axis xy equal tight; colormap parula; colorbar; drawnow;
% 
% %%

clear all

inputs.pixresFP = 9;
inputs.rhoInner = 3;
inputs.rhoOuter = 10;
inputs.angDeg = 160;
inputs.whichSide = 't';
inputs.centering = 'pixel';

[mask1, xis1, etas1] = falco_gen_SW_mask(inputs);
figure(1); imagesc(xis1, etas1, mask1); axis xy equal tight; colormap gray; drawnow;



% Optional inputs
inputs.clockAngDeg = 270;
inputs.FOV = 11;
inputs.whichSide = 'both';
inputs.shape = 'square';

[mask2, xis2, etas2] = falco_gen_SW_mask(inputs);
figure(2); imagesc(xis2, etas2, mask2); axis xy equal tight; colormap gray; drawnow;


% Check D shape and rotation
inputs.clockAngDeg = 10;
inputs.FOV = 11;
inputs.whichSide = 'd';
inputs.shape = 'D';

[mask3, xis3, etas3] = falco_gen_SW_mask(inputs);
figure(3); imagesc(xis3, etas3, mask3); axis xy equal tight; colormap gray; drawnow;


% Check rectangle shape and rotation
inputs.clockAngDeg = -10;
inputs.FOV = 15;
inputs.whichSide = 'ud';
inputs.shape = 'rectangle';

[mask5, xis5, etas5] = falco_gen_SW_mask(inputs);
figure(5); imagesc(xis5, etas5, mask5); axis xy equal tight; colormap gray; drawnow;



% Check shifting the pupil and having a rectangular array. 
inputs.xiFOV = 30;
inputs.etaFOV = 11;
inputs.xiOffset = 10;
inputs.etaOffset = -5;
inputs.angDeg = 180;
inputs.clockAngDeg = 0;
inputs.whichSide = 't';
inputs.shape = 'square';

[mask4, xis4, etas4] = falco_gen_SW_mask(inputs);
figure(4); imagesc(xis4, etas4, mask4); axis xy equal tight; colormap gray; drawnow;


