% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
clear

inputs.Nbeam = 500; % number of points across the pupil diameter

% Optional Inputs
inputs.centering = 'pixel';%'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
inputs.clocking = 10; % CCW rotation [degrees]
inputs.xShear = 0.2; % [pupil diameters]
inputs.yShear = -0.3; % [pupil diameters]

pupil = falco_gen_pupil_Keck(inputs);

figure(1); imagesc(pupil); axis xy equal tight; colorbar; colormap gray;


%% Check if it still runs with all optional inputs left out
clear

inputs.Nbeam = 500; % number of points across the pupil diameter
pupil = falco_gen_pupil_Keck(inputs);

figure(2); imagesc(pupil); axis xy equal tight; colorbar; colormap gray; title('Keck Pupil');


%% Proto-test 1

clear

inputs.Nbeam = 200; % number of points across the pupil diameter

% Optional Inputs
inputs.centering = 'pixel';%'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
inputs.clocking = 65; % CCW rotation [degrees]
inputs.xShear = 0.2; % [pupil diameters]
inputs.yShear = -0.3; % [pupil diameters]

pupil65 = falco_gen_pupil_Keck(inputs);

inputs.clocking = 125; % CCW rotation [degrees]
pupil125 = falco_gen_pupil_Keck(inputs);

sumAbsDiff = sum(sum(abs(pupil125-pupil65)));
sumAbsDiff < 1e-7

figure(3); imagesc(pupil65); axis xy equal tight; colorbar; colormap gray;
figure(4); imagesc(pupil125); axis xy equal tight; colorbar; colormap gray;
figure(5); imagesc(pupil125-pupil65); axis xy equal tight; colorbar; colormap gray;


%% Proto-test 2

clear

inputs.Nbeam = 200; % number of points across the pupil diameter

inputs.centering = 'pixel'; % 'interpixel' or 'pixel'
inputs.clocking = 10; % CCW rotation [degrees]
pupilCentered = falco_gen_pupil_Keck(inputs);

inputs.xShear = 0.5; % [pupil diameters]
inputs.yShear = -0.2; % [pupil diameters]
pupilShifted = falco_gen_pupil_Keck(inputs);
pupilRecentered = circshift(pupilShifted, -inputs.Nbeam*[inputs.yShear, inputs.xShear]);

pupilCentered = pad_crop(pupilCentered, size(pupilShifted));

sumAbsDiff = sum(sum(abs(pupilRecentered-pupilCentered)));
sumAbsDiff < 1e-7

figure(6); imagesc(pupilCentered); axis xy equal tight; colorbar; colormap gray;
figure(7); imagesc(pupilShifted); axis xy equal tight; colorbar; colormap gray;
figure(8); imagesc(pupilRecentered); axis xy equal tight; colorbar; colormap gray;
figure(9); imagesc(pupilRecentered-pupilCentered); axis xy equal tight; colorbar; colormap gray;
