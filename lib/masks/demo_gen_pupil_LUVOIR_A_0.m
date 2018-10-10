% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the LUVOIR Design A telescope pupil in
% Matlab using PROPER (old version)
% Coordinates and dimensions of the struts, primary, secondary, and 
% hex segments are from Matthew Bolcar at NASA GSFC.
%
% The aperture is a six-ring hexagonal array, with the inner most ring removed.  
% There are 120 segments.  Segments are 1.15-m flat-to-flat and segment gaps are 6 mm.  
% The struts are each 125 mm wide.  You can see the strut orientation in the 
% attached figure: they?re all aligned vertically.  
% The bottom two struts are separated by 2.69538 m, center-to-center.
%
% Written by A.J. Riggs on 2017-09-07 to generate the LUVOIR pupil. 
%   Values for the geometry were provided by Matthew Bolcar at NASA GSFC.

close all;
clear;

D_ID = 3.55;
D_OD =  12.644;%12.70;%12.644;
%--15.000 meter tip-to-tip OD of the whole hexagonal array.
% Dfactor = 12.644/(13*1.15 + 12*6e-3)
% Nsp = 1000;
% Nfull = 15.000/D_OD

%--USER INPUTS
NapAcross   = 1000;%250;%324;%1500;%250;                  % number of points across FULL usable pupil
centering = 'interpixel';%'even';%'odd'; % 'pixel'/'odd' or 'interpixel'/'even'

addpath ~/Documents/MATLAB/proper_v3.0.1_matlab_22aug17/

Dap = 13*1.15 + 12*6e-3;
dx = Dap/NapAcross;

NapAcross = 2*ceil(1/2*Dap/dx); % minimum even number of points across to fully contain the actual aperture (if interpixel centered)
if(strcmpi(centering,'pixel'))
    Narray = NapAcross + 2; %--number of points across output array. Requires two more pixels when pixel centered.
else
    Narray = NapAcross; %--number of points across output array. Same size as width when interpixel centered.
end

Darray = Narray*dx;

diam = Darray; %Dap ;%15.05 ;%16; % width of the array (m)
wl_dummy   = 1e-6;               % wavelength (m)
bdf = 1; %--beam diameter factor in output array


width_hex = 1.15;
nrings = 6;
hexrad = 2/sqrt(3)*width_hex/2;
hexsep = width_hex + 6e-3; % 6mm segment gap

switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel','even'}
        cshift = -diam/2/Narray; 
    case {'pixel','odd'}
        cshift = 0;
end
strut_width = 125e-3; % meters

%-------- Generate the input pupil for LUVOIR
bm = prop_begin(Darray, wl_dummy, Narray,'beam_diam_fraction',bdf);

% Subtract the inner ring from all the rings
[bm,ap] = prop_hex_wavefront(bm,nrings,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
[~,ap2] = prop_hex_wavefront(bm,1,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
bm.wf = fftshift(ap-ap2);

%--Add the struts
bm = prop_rectangular_obscuration(bm, strut_width, 8*width_hex, 'XC',cshift, 'YC',cshift + diam/4);
bm = prop_rectangular_obscuration(bm, strut_width, 8*width_hex, 'XC',cshift - 2.69538/2, 'YC',cshift - diam/4);
bm = prop_rectangular_obscuration(bm, strut_width, 8*width_hex, 'XC',cshift + 2.69538/2, 'YC',cshift - diam/4);

pupil = fftshift(bm.wf);
% figure(2); imagesc(pupil); axis xy equal tight; title('Input Pupil','Fontsize',20);
% figure(3); imagesc(pupil-rot90(pupil,2)); axis xy equal tight; title('Input Pupil','Fontsize',20);
% figure(4); imagesc(padOrCropEven(pupil,np/2)); axis xy equal tight; title('Input Pupil','Fontsize',20);
figure(4); imagesc(pupil); axis xy equal tight; title('Input Pupil','Fontsize',20);
% figure(5); imagesc(pupil-rot90(pupil,2)); axis xy equal tight; title('Input Pupil','Fontsize',20);


%%
%%
%%

%-------- Generate the stopped-down input pupil for LUVOIR
bm2 = prop_begin(diam, wl, np);

% Subtract the inner ring from all the rings
[bm2,ap] = prop_hex_wavefront(bm2,nrings,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
[~,ap2] = prop_hex_wavefront(bm2,1,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
bm2.wf = fftshift(ap-ap2);

bm2 = prop_circular_aperture( bm2, D_OD/2,'XC',cshift, 'YC',cshift);
bm2 = prop_circular_obscuration(bm2, D_ID/2,'XC',cshift, 'YC',cshift);

%--Add the struts
bm2 = prop_rectangular_obscuration(bm2, strut_width, 8*width_hex, 'XC',cshift, 'YC',cshift + diam/4);
bm2 = prop_rectangular_obscuration(bm2, strut_width, 8*width_hex, 'XC',cshift - 2.69538/2, 'YC',cshift - diam/4);
bm2 = prop_rectangular_obscuration(bm2, strut_width, 8*width_hex, 'XC',cshift + 2.69538/2, 'YC',cshift - diam/4);

pupil2 = fftshift(bm2.wf);
figure(12); imagesc(padOrCropEven(pupil2,np/2)); axis xy equal tight; title('Input Pupil','Fontsize',20);
% figure(13); imagesc(pupil-rot90(pupil,2)); axis xy equal tight; title('Input Pupil','Fontsize',20);
% figure(14); imagesc(padOrCropEven(pupil,np/2)); axis xy equal tight; title('Input Pupil','Fontsize',20);


tr_ratio = sum(sum(abs(pupil2).^2))/sum(sum(abs(pupil).^2))
%%
figure(13); imagesc(padOrCropEven(pupil-pupil2*.5,np/2)); axis xy equal tight; axis off
% figure(23); imagesc(padOrCropEven(pupil2-pupil*.5,np/2)); axis xy equal tight; axis off

figure(23); imagesc(padOrCropEven(pupil2+pupil,np/2)); axis xy equal tight; axis off
% export_fig('fig_LUVOIR_1D_IDOD_overlap.png','-dpng','-r300','-transparent');

figure(4); imagesc(padOrCropEven(pupil,np/2)); axis xy equal tight; axis off;
% export_fig('fig_LUVOIR_A.png','-dpng','-r300','-transparent');

pupil_LUVOIR = padOrCropEven(pupil,N);