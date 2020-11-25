% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

clear all; 

fnum = 32;
lam = 600e-9; % wavelength [meters]
DspotM = 9*8e-6; % focal spot diameter [meters]

charge = 6;
Nbeam = 256;
Npad = 16*Nbeam;
Nout = 2*Nbeam;
useGPU = false;

%% Input pupil
inputs.Npad = Npad;
inputs.Nbeam = Nbeam; % number of points across the pupil diameter
inputs.OD = 1; % Outer radius (fraction of Nbeam) 
inputs.centering = 'pixel';%'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
Epup1 = falco_gen_pupil_Simple( inputs );

figure(1); imagesc(Epup1); axis xy equal tight; colorbar; drawnow;


%% Propagation with FFTs
vortex = falco_gen_vortex_mask(charge, Npad);

diamSpotLamD = DspotM/(fnum*lam);

inputs.pixresFPM = Npad/Nbeam; %--pixels per lambda/D
inputs.rhoInner = diamSpotLamD/2; % radius of inner FPM amplitude spot (in lambda_c/D)
inputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
inputs.FPMampFac = 0; %0.4; % amplitude transmission of inner FPM spot
inputs.centering = 'pixel';

spot = falco_gen_annular_FPM(inputs);
spot = pad_crop(spot, Npad, 'extrapval', 1); 
figure(2); imagesc(angle(vortex)); axis xy equal tight; colorbar; drawnow;
figure(3); imagesc(spot); axis xy equal tight; colorbar; drawnow;


Efoc1 = fftshift(fft2(fftshift(Epup1)))/Npad;
Ifoc1 = abs(Efoc1).^2;
Ifoc1 = Ifoc1/max(Ifoc1(:));
figure(4); imagesc(log10(Ifoc1));  axis xy equal tight; colorbar;

Epup2fft = fftshift(fft2(fftshift(spot.*vortex.*Efoc1)))/Npad;
Ipup2fft = abs(Epup2fft).^2;

%% Propagation with MFTs with low-res and high-res regions

Epup2 = propcustom_mft_Pup2Vortex2Pup_spot(pad_crop(Epup1, Nout), charge, Nbeam/2, 0.3, 5, useGPU, diamSpotLamD);  %--MFTs
Ipup2 = abs(Epup2).^2;


%% Plots

Ipup2 = pad_crop(Ipup2, Nout);
Ipup2fft = pad_crop(Ipup2fft, Nout);

cs = [-7, 0];
figure(5); imagesc(log10(Ipup2), cs); axis xy equal tight; colorbar;
figure(6); imagesc(log10(Ipup2fft), cs); axis xy equal tight; colorbar;


figure(7); imagesc(angle(Epup2)); axis xy equal tight; colorbar;
figure(8); imagesc(angle(pad_crop(Epup2fft, Nout))); axis xy equal tight; colorbar;