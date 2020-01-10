% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear all;

addpath(genpath('~/Repos/falco-matlab/lib/')) %--Add FALCO library to MATLAB path
addpath(genpath('~/Documents/MATLAB/PROPER/')) %--Add PROPER library to MATLAB path

%--Load true pupil
pupil0 = imread('2019-10-09 CGI entrance pupil 8k binary.png');
pupilFromFile = double(pupil0)/double(max(pupil0(:)));
Narray = size(pupilFromFile,1);

%--Generate pupil representation in FALCO
Nbeam = 2*4022.0;
centering = 'pixel';
changes.dummy = 1;
changes.xShear = -0.5/Nbeam; 
changes.yShear = -2/Nbeam;
pupilFromFALCO = falco_gen_pupil_WFIRST_CGI_20191009(Nbeam, centering, changes);
pupilFromFALCO = padOrCropEven(pupilFromFALCO, Narray);

figure(1); imagesc(pupilFromFile); axis xy equal tight; colorbar; drawnow;
figure(2); imagesc(pupilFromFALCO); axis xy equal tight; colorbar; drawnow;
figure(3); imagesc(pupilFromFile - pupilFromFALCO, [-1 1]); axis xy equal tight; colorbar; colormap gray; drawnow;
