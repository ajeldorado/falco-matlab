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
% Written on 2018-05-16 by A.J. Riggs to generate the HabEx B pupil. 

% function pupil = falco_gen_pupil_HabEx_B_PROPER(inputs,varargin)

%-------------------
%--FOR DEBUGGING ONLY
clear;
addpath ~/Repos/FALCO/lib/PROPER/
addpath ~/Repos/FALCO/lib/utils/
inputs.Nbeam = 1000;
inputs.gap_width = 25e-3; % (meters)
inputs.centering = 'pixel';
%-------------------

%--USER INPUTS
Nbeam   = inputs.Nbeam; % number of points across FULL usable pupil
centering = inputs.centering;% 'pixel' or 'interpixel' centering of the array
gap_width = inputs.gap_width;


%--Do not change
Dap = 4.0; % Circumscribing aperture diameter (meters)
dx = Dap/Nbeam;
hexradius = .2374*4; %--Radius of circumscribing circle for the inner hex (NOT INCLUDING SEGMENT GAP) (meters)

if(strcmpi(centering,'pixel'))
    Narray = Nbeam + 2; %--number of points across output array. Requires two more pixels when pixel centered.
else
    Narray = Nbeam;
end

%--Values for PROPER
Darray = Narray*dx; % width of the array (m)
wl_dummy   = 1e-6;               % wavelength (m)
bdf = 1; %--beam diameter factor in output array

%--Centering of the aperture on the array
switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel'}
        cshift = -Darray/2/Narray; 
    case {'pixel'}
        cshift = 0;
    otherwise
        error('falco_gen_pupil_LUVOIR_A_0.m: Error! Centering must be either pixel or interpixel.')
end

%-------- Generate the input pupil for LUVOIR with PROPER
bm = prop_begin(Darray, wl_dummy, Narray,'beam_diam_fraction',bdf);


%--OUTER CIRCLE
ra_OD = Dap/2;
cx_OD = 0 + cshift;
cy_OD = 0 + cshift;
bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);

%--Inner Hex Ring
hexOut = prop_polygon( bm, 6, hexradius+gap_width , 'XC', cshift  , 'YC', cshift, 'DARK' , 'ROTATION', 30 ); 
hexIn = prop_polygon( bm, 6, hexradius , 'XC', cshift  , 'YC', cshift , 'ROTATION', 30 ); 
bm.wf = bm.wf.*fftshift(hexOut+hexIn);


%--Rectangular Gaps:
buffer = 2*gap_width;
length_rect = (Dap - 2*hexradius)/2 + buffer;
R_rect_cent = length_rect/2 + hexradius ;

%--12 O'clock Gap
Nrect = 6;
for ii=1:Nrect
    angDeg = (ii-1)*60;
    bm = prop_rectangular_obscuration(bm, gap_width, length_rect, 'XC',cshift + R_rect_cent*sind(angDeg), 'YC',cshift + R_rect_cent*cosd(angDeg),'ROTATION',-angDeg );
end

pupil = fftshift(bm.wf);

% figure(1); imagesc(hexOut); axis xy equal tight; title('Hexagon','Fontsize',20); colorbar;
% figure(2); imagesc(pupil.*(hexOut+hexIn)); axis xy equal tight; title('Pupil','Fontsize',20); colorbar;
% figure(3); imagesc(pupil); axis xy equal tight; title('Pupil','Fontsize',20); colorbar;
figure(3); imagesc(pupil); axis xy equal tight; axis off; colormap gray;
% export_fig('fig_HabExB_pupil.png','-dpng','-r200')

%%
% %--Check centering for case 'interpixel'
% figure(5); imagesc(pupil-fliplr(pupil)); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;
% figure(6); imagesc(pupil-rot90(pupil,2)); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;
% 
% %--Check centering for case 'pixel'
% pupilCrop = pupil(2:end,2:end);
% figure(15); imagesc(pupilCrop-fliplr(pupilCrop)); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;
% figure(16); imagesc(pupilCrop-rot90(pupilCrop,2)); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;

%%

% %% Compare with the actual Aperture
% addpath ~/Repos/FALCO/lib/masks/HabExB/
% 
% fn = 'OAP1_6mmGap_5mmRadii.png';
% 
% temp = imread(fn);
% mask = temp(:,:,1);
% mask = 1-mask;
% sum0 = sum(mask(:))
% % figure(10); imagesc(mask); axis xy equal tight; title('Pupil 0','Fontsize',20); colorbar;
% 
% mask = mask(12:end,:);
% sum1 = sum(mask(:))
% % figure(11); imagesc(mask); axis xy equal tight; title('Pupil 0','Fontsize',20); colorbar;
% 
% mask = mask(1:3817,:);
% sum2 = sum(mask(:))
% % figure(12); imagesc(mask); axis xy equal tight; title('Pupil 0','Fontsize',20); colorbar;
% 
% mask = mask(:,66:end);
% sum3 = sum(mask(:))
% % figure(13); imagesc(mask); axis xy equal tight; title('Pupil 0','Fontsize',20); colorbar;
% 
% mask = mask(:,1:3818);
% sum4 = sum(mask(:))
% % figure(14); imagesc(mask); axis xy equal tight; title('Cropped Pupil','Fontsize',20); colorbar;
% 
% Nx = size(mask,2);
% Ny = size(mask,1);
% 
% xs = ( -(Nx-1)/2:(Nx-1)/2 )/Nx;
% ys = ( -(Ny-1)/2:(Ny-1)/2 )/Ny;
% 
% figure(15); imagesc(xs,ys,mask); axis xy equal tight; title('Cropped Pupil','Fontsize',20); colorbar;
% 
% %%


% end %---END OF FUNCTION

