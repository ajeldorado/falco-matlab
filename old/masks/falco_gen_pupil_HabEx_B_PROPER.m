% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function pupil = falco_gen_pupil_LUVOIR_A_0(inputs,varargin)
%
%--Function to generate an approximate HabEx B segmented aperture. The
%  rounded corners at vertices are left out because they are harder to
%  generate. The segment gap size can be specified.
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Written on 2018-05-16 by A.J. Riggs to generate the HabEx B pupil. 

function pupil = falco_gen_pupil_HabEx_B_PROPER(inputs,varargin)

% %-------------------
% %--FOR DEBUGGING ONLY
% clear;
% addpath ~/Repos/FALCO/lib/PROPER/
% addpath ~/Repos/FALCO/lib/utils/
% inputs.Nbeam = 1000;
% inputs.wGap = 6e-3; % (meters)
% inputs.centering = 'pixel';
% %-------------------

%--USER INPUTS
Nbeam   = inputs.Nbeam; % number of points across FULL usable pupil
centering = inputs.centering;% 'pixel' or 'interpixel' centering of the array
wGap = inputs.wGap;

%--Do not change
Dmask = 4.0; % Circumscribing aperture diameter (meters)
dx = Dmask/Nbeam;
hexradius = .2374*4; %--Radius of circumscribing circle for the inner hex (NOT INCLUDING SEGMENT GAP) (meters)

if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam+1); %--number of points across output array. Requires two more pixels when pixel centered.
else
    Narray = ceil_even(Nbeam);
end

%--Values for PROPER
Darray = Narray*dx; % width of the array (m)
wl_dummy   = 1e-6;               % wavelength (m)
bdf = Nbeam/Narray; %--beam diameter factor in output array

%--Centering of the aperture on the array
switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
    otherwise
        error('falco_gen_pupil_HabEx_B_PROPER.m: Error! Centering must be either pixel or interpixel.')
end

%-------- Generate the input pupil for LUVOIR with PROPER
bm = prop_begin(Dmask, wl_dummy, Narray,'beam_diam_fraction',bdf);

%--OUTER CIRCLE
ra_OD = Dmask/2;
cx_OD = 0 + cshift;
cy_OD = 0 + cshift;
bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);

%--Inner Hex Ring
hexOut = prop_polygon( bm, 6, hexradius+wGap , 'XC', cshift  , 'YC', cshift, 'DARK' , 'ROTATION', 30 ); 
hexIn = prop_polygon( bm, 6, hexradius , 'XC', cshift  , 'YC', cshift , 'ROTATION', 30 ); 
bm.wf = bm.wf.*fftshift(hexOut+hexIn);

%--Rectangular Gaps:
buffer = 2*wGap;
length_rect = (Dmask - 2*hexradius)/2 + buffer;
R_rect_cent = length_rect/2 + hexradius ;

%--12 O'clock Gap
Nrect = 6;
for ii=1:Nrect
    angDeg = (ii-1)*60;
    bm = prop_rectangular_obscuration(bm, wGap, length_rect, 'XC',cshift + R_rect_cent*sind(angDeg), 'YC',cshift + R_rect_cent*cosd(angDeg),'ROTATION',-angDeg );
end

pupil = fftshift(bm.wf);

end %---END OF FUNCTION