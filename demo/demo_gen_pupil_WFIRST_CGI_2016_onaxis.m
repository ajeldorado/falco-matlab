% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear all;

%%
clear all
% Test magnification, rotation, and shear of pupil
DeltaY = 2; % pixels
Nbeam = 1000;
centering = 'interpixel';
changes.magFac = 0.7;
changes.yShear = -DeltaY/Nbeam;

changes.clock_deg = 0;
pupilA = falco_gen_pupil_WFIRST_2016_onaxis(Nbeam, centering, changes);
pupilBprime = circshift(rot90(pupilA, 2), [-2*DeltaY, 0]);


changes.clock_deg = 180;
pupilB = falco_gen_pupil_WFIRST_2016_onaxis(Nbeam, centering, changes);

figure(1); imagesc(pupilB); axis xy equal tight; colorbar; colormap gray; drawnow;
figure(2); imagesc(pupilBprime); axis xy equal tight; colorbar; colormap gray; drawnow;
figure(3); imagesc(pupilB - pupilBprime); axis xy equal tight; colorbar; colormap gray; drawnow;


% figure(1); imagesc(pupilA - fliplr(pupilA)); axis xy equal tight; colorbar; colormap gray; drawnow;


%% Lyot stop generation

Nbeam = 309;
centering = 'pixel';

clear changes
changes.flagLyot = true;
changes.ID = 0.50;
changes.OD = 0.80;
changes.wStrut = 0.036;
lyot = falco_gen_pupil_WFIRST_2016_onaxis(Nbeam, centering, changes);
lyotCropped = lyot(2:end, 2:end);

figure(4); imagesc(lyot); axis xy equal tight; colorbar; colormap gray; drawnow;
figure(5); imagesc(lyotCropped - fliplr(lyotCropped)); axis xy equal tight; colorbar; colormap gray; drawnow;


%% Test 90 degree pupil rotation
clear changes 

centering = 'interpixel';
changes.xShear = 0;
changes.yShear = 0;
changes.clock_deg = 90;
pupilB = falco_gen_pupil_WFIRST_2016_onaxis(Nbeam, centering, changes);

changes.clock_deg = 0;
pupilA = falco_gen_pupil_WFIRST_2016_onaxis(Nbeam, centering, changes);

figure(11); imagesc(pupilB); axis xy equal tight; colorbar; colormap gray; drawnow;

figure(12); imagesc(pupilB - rot90(pupilA, -1)); axis xy equal tight; colorbar; colormap gray; drawnow;
