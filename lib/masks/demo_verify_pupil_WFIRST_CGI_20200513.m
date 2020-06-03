% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear all;

% addpath(genpath('~/Repos/falco-matlab/lib/')) %--Add FALCO library to MATLAB path
% addpath(genpath('~/Documents/MATLAB/PROPER/')) %--Add PROPER library to MATLAB path

%--Load true pupil
pupil0 = imread('pupil_CGI-20200513_8k_binary_noF.png');
% pupil0 = rot90(pupil0, 2);
pupilFromFile = double(pupil0)/double(max(pupil0(:)));
Narray = size(pupilFromFile,1);

%--Generate pupil representation in FALCO
Nbeam = 2*4027.25;
centering = 'pixel';
changes.dummy = 1;
changes.xShear = -0.5/Nbeam; 
changes.yShear = -52.85/Nbeam;
pupilFromFALCO = falco_gen_pupil_WFIRST_CGI_20200513(Nbeam, centering, changes);
pupilFromFALCO = pad_crop(pupilFromFALCO, Narray);

diff = pupilFromFile - pupilFromFALCO;
errorSum = sum(sum(abs(diff)));
fprintf('sum of mismatch = %.2f\n', errorSum)
fprintf('sum of underpadded pixels = %.2f\n', sum(diff(diff<-0.5)))
fprintf('sum of overpadded pixels  = %.2f\n', sum(diff(diff>0.5)))


figure(1); imagesc(pupilFromFile); axis xy equal tight; colorbar; drawnow;
figure(2); imagesc(pupilFromFALCO); axis xy equal tight; colorbar; drawnow;
figure(3); imagesc(pupilFromFile - pupilFromFALCO, [-1 1]); axis xy equal tight; colorbar; colormap gray; drawnow;
