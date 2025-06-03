% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%
% Script to demonstrate applying aberrations at a focal plane with
% propcustom_mft_apply_focal_errors_babinet()
% -------------------------------------------------------------------------
clear; 

inputs.Npad = 700;
inputs.Nbeam = 500; % number of points across the pupil diameter
inputs.OD = 1; % Outer radius (fraction of Nbeam) 
inputs.ID = 0.2;% Inner radius (zero if you want an off-axis telescope)
inputs.angStrut = [0 120 240];%Angles of the struts (deg)
inputs.wStrut = 0.01; % Width of the struts (fraction of pupil diam.)
inputs.centering = 'pixel';%'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
inputs.stretch = 1;% Stretch the horizontal axis to create elliptical beam 

% Optional Inputs
inputs.clocking = 15; % CCW rotation. Doesn't work with flag HG. [degrees]
inputs.xShear = 0.1;

pupil = falco_gen_pupil_Simple( inputs );

figure(1);
imagesc(pupil);
axis xy equal tight; colorbar;
colormap gray;


Nfoc = 200;
resErrorMap = 4; % pixels per lambda/D
zeroMeanNormNoise = (rand(Nfoc, Nfoc) - 0.5)/0.2887;
% EerrorMap = exp(1j*zeroMeanNormNoise*0.8/360);
EerrorMap = 1 + 0.1*zeroMeanNormNoise;
beamDiamPix = inputs.Nbeam;


%% Try making a Lyot spot as the error pattern
% Parameters
resErrorMap = 6;
beamDiamPix = inputs.Nbeam;
r = 10; %2.8;  % pixels per lambda/D
% Create array with central circular region
[x, y] = meshgrid(1:beamDiamPix, 1:beamDiamPix);
EerrorMap = sqrt((x - ceil(beamDiamPix/2)).^2 + (y - ceil(beamDiamPix/2)).^2) > r*resErrorMap;


figure(11);
imagesc(abs(EerrorMap));
axis xy equal tight; colorbar;
colormap gray;

figure(12);
imagesc(angle(EerrorMap));
axis xy equal tight; colorbar;
colormap hsv;


% The direct approach has vignetting issues without using large arrays:
EpupNew = propcustom_mft_apply_focal_errors(pupil, EerrorMap, resErrorMap, beamDiamPix);

figure(2);
imagesc(abs(EpupNew));
axis xy equal tight; colorbar;
colormap gray;


% Babinet's principle allows for a smaller array of errors at the focal
% plane:
EpupNewB = propcustom_mft_apply_focal_errors_babinet(pupil, EerrorMap, resErrorMap, beamDiamPix);

figure(3);
imagesc(abs(EpupNewB));
axis xy equal tight; colorbar;
colormap gray;

figure(4);
imagesc(angle(EpupNewB));
axis xy equal tight; colorbar;
colormap hsv;

figure(5);
imagesc(abs(EpupNewB) - pupil);
axis xy equal tight; colorbar;
colormap gray;
