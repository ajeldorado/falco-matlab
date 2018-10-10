% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the LUVOIR Design A, Aperture 5, telescope pupil in
%   Matlab using PROPER
% Coordinates and dimensions of the primary, secondary, and hex segments
%   are from Matthew Bolcar (NASA GSFC).
% Coordinates and dimenstions of the secondary mirror support struts were a
%   best-fit match by A.J. Riggs by matching PROPER-made rectangles to the 
%   pupil file from Matthew Bolcar (NASA GSFC).
%
% Modified on 2018-02-25 by A.J. Riggs to be for LUVOIR A aperture 5. 
% Written on  2017-09-07 by A.J. Riggs to generate the first proposed LUVOIR pupil. 
%   Values for the geometry were provided by Matthew Bolcar at NASA GSFC.
%
%--Coordinates of hex segments to skip:
% 1 13 114 115 126 127
% 1 12 113 114 125 126

% %--FOR DEBUGGING ONLY
clear;
addpath ../PROPER/
inputs.Nbeam = 1000;
inputs.magfacD = 1;
inputs.strut_width = 1/100;


mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);

pupil = mask;%fftshift(bm.wf);
figure(2); imagesc(pupil); axis xy equal tight; title('Input Pupil','Fontsize',20);

% 
% % function mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs,varargin)
% 
% %--FOR DEBUGGING ONLY
% clear;close all;
% % cd /Users/ajriggs/Repos/FALCO/DoNotPush
% addpath ../PROPER/
% inputs.Nbeam = 1000;
% inputs.magfacD = 1;
% inputs.centering = 'interpixel';
% flagRot180deg = false;
% inputs.strutWidth = 0.2/100;
% % % inputs.Dbeam = 14.9760; % (m)
% % % inputs.xshift = 0;
% % % inputs.yshift = 0;
% %-------------------
% 
% % %--Set default values of input parameters
% % flagRot180deg = false;
% % %--Look for Optional Keywords
% % icav = 0;                     % index in cell array varargin
% % while icav < size(varargin, 2)
% %     icav = icav + 1;
% %     switch lower(varargin{icav})
% %         case {'rot180'}
% %             flagRot180deg = true; % For even arrays, beam center is in between pixels.
% %         otherwise
% %             error('falco_gen_pupil_WFIRSTcycle6_mag_rot_trans: Unknown keyword: %s\n', ...
% %             varargin{icav});
% %     end
% % end
% 
% %--Centering of array: 'pixel' or 'interpixel'
% if(isfield(inputs,'centering'))
%     centering = inputs.centering;
% else %--Default to pixel centering
%     centering = 'pixel';
% end
% 
% %--Magnification factor of the pupil diameter
% if(isfield(inputs,'magfacD')) 
%     magfacD = inputs.magfacD;
% else
%     magfacD = 1;
% end
% 
% Nbeam   = inputs.Nbeam;     % number of points across the incoming beam  
% % magfacD = inputs.magfacD; %magnification factor of the pupil diameter
% % Dbeam = inputs.Dbeam; %--diameter of the beam at the mask (meters)
% % % Narray = inputs.Narray;   % number of points across output array
% % clock_deg = inputs.clock_deg; % clocking angle of the pupil (in degrees)
% % xshift = inputs.xshift; % translation in x of pupil (in diameters)
% % yshift = inputs.yshift; % translation in y of pupil (in diameters)
% 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% % %
% 
% 
% 
% 
% %--USER INPUTS
% Nap   = Nbeam;%250;%324;%1500;%250;                  % number of points across FULL usable pupil
% % magfac_vec = 1 + (-1:0.1:1)/100;
% % magfac_vec = 1 + (-2:0.1:2)/100;
% 
% 
% % Nmag = 1; %length(magfac_vec);
% 
% 
% width_hex0 = 1.242; %-- flat-to-flat (m)
% hexgap0 = 6e-3; % (m)
% Dap = (12*width_hex0 + 12*hexgap0);
% dx = Dap/Nap;
% dx_drawing = 1.242/158; % (m) %--In actual drawing, 158 pixels across the 1.242m, so dx_pixel =
% %dx_drawing. strut in drawing is 19 pixels across, so dx_strut =
% %19*dx_drawing = 0.1494m
% 
% Narray = 2*ceil(1/2*Dap/dx/magfacD); % minimum even number of points across to fully contain the actual aperture (if interpixel centered)
% % Narray = 2*ceil(1/2*Dap/dx/max(magfac_vec)); % minimum even number of points across to fully contain the actual aperture (if interpixel centered)
% if(strcmpi(centering,'pixel') || strcmpi(centering,'odd'))
%     Narray = Narray + 2; %--number of points across output array. Requires two more pixels when pixel centered.
% else
% %     Narray = Narray; %--number of points across output array. Same size as width when interpixel centered.
% end
% Darray = Narray*dx;
% % pupil_sum = zeros(Narray,Narray);
% 
% %--For PROPER 
% diam = Darray; %Dap;%15.05 ;%16; % width of the array (m)
% wl_dummy   = 1e-6;               % wavelength (m)
% bdf = 1; %--beam diameter factor in output array
% 
% 
%         
% dx_t = 0;%xyvec(ti,1);
% dy_t = 0;%xyvec(ti,2);
% 
% magfac = magfacD; %magfac_vec(mi);
% 
% width_hex = magfac*width_hex0; %-- flat-to-flat (m)
% nrings = 6;
% hexrad = 2/sqrt(3)*width_hex/2;
% hexgap = magfac*hexgap0; % (m)
% hexsep = width_hex + hexgap; % distance from center to center of neighboring segments
% 
% 
% 
% 
% switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
%     case {'interpixel','even'}
%         cshift = -diam/2/Narray; 
%     case {'pixel','odd'}
%         cshift = 0;
% end
% 
% strut_width0 = 19*dx_drawing*magfac; %%%%125e-3; % meters
% strut_width = strut_width*diam; %--now in meters
% 
% 
% %-------- Generate the input pupil for LUVOIR
% bm = prop_begin(Darray, wl_dummy, Narray,'beam_diam_fraction',bdf);
% 
% % Subtract the inner ring from all the rings
% [ap] = falco_hex_aperture_LUVOIR_A_5(bm,nrings,hexrad,hexsep,'XC',cshift-dx_t,'YC',cshift-dy_t,'DARKCENTER'); %--Official Matlab PROPER from August 2017
% % [~,ap2] = prop_hex_wavefront(bm,0,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
% % bm.wf = fftshift(ap-ap2);
% 
% % %--Add the struts
% bm = prop_rectangular_obscuration(bm, strut_width, 7*width_hex, 'XC',cshift-dx_t, 'YC',cshift-dy_t + magfac*diam/4);
% 
% 
% len_1a = 2*width_hex - 12*dx_drawing;
% bm = prop_rectangular_obscuration(bm, strut_width, len_1a, 'XC',cshift-dx_t + (hexrad-0.5*strut_width0), 'YC',cshift-dy_t - len_1a/2.);
% bm = prop_rectangular_obscuration(bm, strut_width, len_1a, 'XC',cshift-dx_t - (hexrad-0.5*strut_width0), 'YC',cshift-dy_t - len_1a/2.);
% 
% len_1b = 3.75*width_hex;
% bm = prop_rectangular_obscuration(bm, strut_width, len_1b, 'XC',cshift-dx_t + 1.25*hexrad*2, 'YC',cshift-dy_t - 3.5*width_hex,'ROT',30);
% bm = prop_rectangular_obscuration(bm, strut_width, len_1b, 'XC',cshift-dx_t - 1.25*hexrad*2, 'YC',cshift-dy_t - 3.5*width_hex,'ROT',-30);
% 
% % bm = prop_rectangular_obscuration(bm, strut_width, 2*width_hex, 'XC',cshift, 'YC',cshift - diam/8);
% 
% % bm = prop_rectangular_obscuration(bm, strut_width, 8*width_hex, 'XC',cshift - 2.69538/2, 'YC',cshift - diam/4);
% % bm = prop_rectangular_obscuration(bm, strut_width, 8*width_hex, 'XC',cshift + 2.69538/2, 'YC',cshift - diam/4);
% mask = ifftshift(abs(bm.wf)).*ap;
% 
% figure(2); imagesc(mask); axis xy equal tight; title('Input Pupil','Fontsize',20);
%     
% % switch centering
% %     case{'interpixel'}
% %         figure(4); imagesc(pupil-fliplr(pupil)); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;
% %     case{'pixel'}      
% %         figure(5); imagesc(pupil(2:end,2:end)-fliplr(pupil(2:end,2:end))); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;
% % end
% 
% 
% 
% 
% %%
% %%
% %%
% 
% D_OD = Darray*13.7/15.2;
% D_ID = (2*0.05615)*Darray*13.7/15.2;
% 
% 
% % % % Darray = Narray*dxiUL; %--width of array in lambda_c/D
% % % % diam = Darray;  
% wl_dummy   = 1e-6;              % wavelength (m); Dummy value--no propagation here, so not used.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% %--INITIALIZE PROPER
% bdf = 1; %--beam diameter factor in output array
% bm2 = prop_begin(diam, wl_dummy, Narray,'beam_diam_fraction',bdf);
% 
% 
% 
% 
% %-------- Generate the stopped-down input pupil for LUVOIR
% % bm2 = prop_begin(diam, wl, np);
% 
% % % Subtract the inner ring from all the rings
% % [bm2,ap] = prop_hex_wavefront(bm2,nrings,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
% % [~,ap2] = prop_hex_wavefront(bm2,1,hexrad,hexsep,'XCENTER',cshift,'YCENTER',cshift); %--Official Matlab PROPER from August 2017
% % bm2.wf = fftshift(ap-ap2);
% 
% bm2 = prop_circular_aperture( bm2, D_OD/2,'XC',cshift, 'YC',cshift);
% bm2 = prop_circular_obscuration(bm2, D_ID/2,'XC',cshift, 'YC',cshift);
% 
% % %--Add the struts
% % bm2 = prop_rectangular_obscuration(bm2, strut_width, 8*width_hex, 'XC',cshift, 'YC',cshift + diam/4);
% % bm2 = prop_rectangular_obscuration(bm2, strut_width, 8*width_hex, 'XC',cshift - 2.69538/2, 'YC',cshift - diam/4);
% % bm2 = prop_rectangular_obscuration(bm2, strut_width, 8*width_hex, 'XC',cshift + 2.69538/2, 'YC',cshift - diam/4);
% 
% LS = fftshift(bm2.wf);
% figure(12); imagesc(LS); axis xy equal tight; title('Input Pupil','Fontsize',20);
% figure(13); imagesc(mask+LS); axis xy equal tight; title('Input Pupil','Fontsize',20);
% % figure(13); imagesc(pupil-rot90(pupil,2)); axis xy equal tight; title('Input Pupil','Fontsize',20);
% % figure(14); imagesc(padOrCropEven(pupil,np/2)); axis xy equal tight; title('Input Pupil','Fontsize',20);
% 
% 
% tr_ratio = sum(sum(abs(LS).^2))/sum(sum(abs(pupil).^2))
% %%
% figure(13); imagesc(padOrCropEven(pupil-LS*.5,np/2)); axis xy equal tight; axis off
% % figure(23); imagesc(padOrCropEven(pupil2-pupil*.5,np/2)); axis xy equal tight; axis off
% 
% figure(23); imagesc(padOrCropEven(LS+pupil,np/2)); axis xy equal tight; axis off
% % export_fig('fig_LUVOIR_1D_IDOD_overlap.png','-dpng','-r300','-transparent');
% 
% figure(4); imagesc(padOrCropEven(pupil,np/2)); axis xy equal tight; axis off;
% % export_fig('fig_LUVOIR_A.png','-dpng','-r300','-transparent');
% 
% pupil_LUVOIR = padOrCropEven(pupil,N);
% 
% 
% 
% % end
% 
% 
