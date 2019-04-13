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
%  -FOV: minimum desired field of view (in lambda_c/D)
%  -shape: 'square' makes a square. Omitting makes a circle. 
%  -clockAngDeg: Dark hole rotation about the z-axis (deg)
%
%--OUTPUTS:
% maskSW: rectangular, even-sized, binary-valued software mask
% xis: vector of coordinates along the horizontal axis (in lambda_c/D)
% etas: : vector of coordinates along the vertical axis (in lambda_c/D)

clear; 
addpath('../utils/');

input.pixresFP = 9;
input.rhoInner = 3;
input.rhoOuter = 10;
input.angDeg = 160;
input.whichSide = 'r';
input.clockAngDeg = 0;
input.centering = 'pixel';
input.FOV = 512/input.pixresFP/2;
% input.shape = 'square';
mask = falco_gen_SW_mask(input);

%%

figure(1);
imagesc(mask);
axis image;
set(gca,'ydir','normal');
