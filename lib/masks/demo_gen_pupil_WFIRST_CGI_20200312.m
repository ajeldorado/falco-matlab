% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear all;

% addpath(genpath('~/Repos/falco-matlab/lib/')) %--Add FALCO library to MATLAB path
% addpath(genpath('~/Documents/MATLAB/PROPER/')) %--Add PROPER library to MATLAB path

Nbeam = 1045;
centering = 'pixel';
changes.clock_deg = 15;
changes.xShear = 0.1;
changes.yShear = 0.2;
changes.magFac = 0.5;
pupil = falco_gen_pupil_WFIRST_CGI_20200312(Nbeam, centering, changes);

figure(1); imagesc(pupil); axis xy equal tight; colorbar; drawnow;

%% Lyot stop mode

Nbeam = 309;
centering = 'pixel';

clear changes
changes.flagLyot = true;
changes.ID = 0.50;
changes.OD = 0.80;
LS = falco_gen_pupil_WFIRST_CGI_20200312(Nbeam, centering, changes);

figure(2); imagesc(LS); axis xy equal tight; colorbar; drawnow;

croppedLS = LS(2:end,2:end);
figure(3); imagesc(croppedLS - fliplr(croppedLS)); axis xy equal tight; colorbar; drawnow;
