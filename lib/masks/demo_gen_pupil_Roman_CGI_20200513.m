% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear all;

DeltaY = 2; % pixels
Nbeam = 1000;
centering = 'interpixel';
changes.magFac = 0.7;
changes.yShear = -DeltaY/Nbeam;

changes.clock_deg = 0;
pupilA = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);


changes.clock_deg = 180;
pupilC = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);

figure(1); imagesc(pupilA); axis xy equal tight; colorbar; colormap gray; drawnow;

figure(5); imagesc(pupilC - circshift(rot90(pupilA, 2), [-2*DeltaY, 0])); axis xy equal tight; colorbar; colormap gray; drawnow;

figure(3); imagesc(pupilC); axis xy equal tight; colorbar; colormap gray; drawnow;
figure(6); imagesc(rot90(pupilA, 2)); axis xy equal tight; colorbar; colormap gray; drawnow;


changes.xShear = 0;
changes.yShear = 0;
changes.clock_deg = 90;
pupilB = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);

changes.clock_deg = 0;
pupilA = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);

figure(2); imagesc(pupilB); axis xy equal tight; colorbar; colormap gray; drawnow;

figure(4); imagesc(pupilB - rot90(pupilA, -1)); axis xy equal tight; colorbar; colormap gray; drawnow;
