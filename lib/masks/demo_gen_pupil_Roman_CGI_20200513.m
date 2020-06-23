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


%%
Nbeam = 1045;
centering = 'pixel';
changes.clock_deg = 15;
changes.xShear = 0.1;
changes.yShear = 0.2;
changes.magFac = 0.5;
pupil = falco_gen_pupil_WFIRST_CGI_20200513(Nbeam, centering, changes);

pupilB = falco_gen_pupil_WFIRST_CGI_20200513_cleaner(Nbeam, centering, changes);

figure(1); imagesc(pupil); axis xy equal tight; colorbar; colormap gray; drawnow;
figure(2); imagesc(pupilB); axis xy equal tight; colorbar; colormap gray; drawnow;

diff = pupil-pupilB;
% diff(diff>0) = 1;
% diff(diff<0) = -1;

figure(3); imagesc(diff); axis xy equal tight; colormap jet; colorbar; drawnow;
fprintf('Sum of absolute difference of pupils = %.3f\n', sum(sum(abs(pupil-pupilB))));

% %% Lyot stop mode

% Nbeam = 309;
% centering = 'pixel';
% 
% clear changes
% changes.flagLyot = true;
% changes.ID = 0.50;
% changes.OD = 0.80;
% LS = falco_gen_pupil_WFIRST_CGI_20200513(Nbeam, centering, changes);
% 
% figure(2); imagesc(LS); axis xy equal tight; colorbar; drawnow;
% 
% croppedLS = LS(2:end,2:end);
% figure(3); imagesc(croppedLS - fliplr(croppedLS)); axis xy equal tight; colorbar; drawnow;
