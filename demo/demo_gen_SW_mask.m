% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to demonstrate the usage cases of falco_gen_SW_mask.
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

inputs.pixresFP = 9;
inputs.rhoInner = 3;
inputs.rhoOuter = 10;
inputs.angDeg = 160;
inputs.whichSide = 't';
inputs.centering = 'pixel';

[mask1, xis1, etas1] = falco_gen_SW_mask(inputs);
figure(1); imagesc(xis1, etas1, mask1); axis xy equal tight; colormap gray; drawnow;


% Optional inputs
inputs.FOV = 11;
inputs.whichSide = 'ud';
inputs.shape = 'square';

[mask2, xis2, etas2] = falco_gen_SW_mask(inputs);
figure(2); imagesc(xis2, etas2, mask2); axis xy equal tight; colormap gray; drawnow;


% Check D shape and rotation
inputs.clockAngDeg = 0;
inputs.FOV = 11;
inputs.whichSide = 'top';
inputs.shape = 'D';

[mask3, xis3, etas3] = falco_gen_SW_mask(inputs);
figure(3); imagesc(xis3, etas3, mask3); axis xy equal tight; colormap gray; drawnow;


% Check rectangle shape and rotation
inputs.clockAngDeg = -10;
inputs.FOV = 15;
inputs.whichSide = 'ud';
inputs.shape = 'rectangle';

[mask4, xis4, etas4] = falco_gen_SW_mask(inputs);
figure(4); imagesc(xis4, etas4, mask4); axis xy equal tight; colormap gray; drawnow;


% Check shifting the mask and having a rectangular array. 
inputs.xiFOV = 30;
inputs.etaFOV = 11;
inputs.xiOffset = 10;
inputs.etaOffset = -5;
inputs.angDeg = 180;
inputs.clockAngDeg = 0;
inputs.whichSide = 't';
inputs.shape = 'square';

[mask5, xis5, etas5] = falco_gen_SW_mask(inputs);
figure(5); imagesc(xis5, etas5, mask5); axis xy equal tight; colormap gray; drawnow;


% Check making a square offset from the star 
clear inputs
inputs.pixresFP = 9;
inputs.rhoInner = 0;
inputs.rhoOuter = 2.5;
inputs.xiFOV = 16;
inputs.etaFOV = 20;
inputs.xiOffset = 4;
inputs.etaOffset = 0;
inputs.angDeg = 180;
inputs.clockAngDeg = 0;
inputs.whichSide = 'lr';
inputs.shape = 'square';

[mask6, xis6, etas6] = falco_gen_SW_mask(inputs);
figure(6); imagesc(xis6, etas6, mask6); axis xy equal tight; colormap gray; drawnow;


clear;

inputs.pixresFP = 9;
inputs.rhoInner = 3;
inputs.rhoOuter = 10;
inputs.angDeg = 180;
inputs.clockAngDeg = 270;
inputs.whichSide = 't';
inputs.centering = 'pixel';

[mask11, xis11, etas11] = falco_gen_SW_mask(inputs);
figure(11); imagesc(xis11, etas11, mask11); axis xy equal tight; colormap gray; drawnow;

inputs.radius_erode = 1;
[mask12, xis12, etas12] = falco_gen_SW_mask(inputs);
figure(12); imagesc(xis11, etas11, mask11+mask12); axis xy equal tight; colormap parula; drawnow;

