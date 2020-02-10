% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

clear all; 

inputs.Nbeam = 200; % max aperture diameter in samples
inputs.Narray = 2^ceil(log2(inputs.Nbeam+2)); % Number of samples across in square output array

inputs.primaryRadiusX = 0.3; % x-radius of ellipse [pupil diameters]
inputs.primaryRadiusY = 0.5; % y-radius of ellipse [pupil diameters]
inputs.primaryClockingDegrees = 30; % clocking of the primary [degrees]

inputs.secondaryRadiusX = 0.2; % x-radius of ellipse [pupil diameters]
inputs.secondaryRadiusY = 0.1; % y-radius of ellipse [pupil diameters]
inputs.secondaryClockingDegrees = -30; % clocking of the secondary [degrees]

inputs.wStrut = 0.02;
inputs.angStrut = [-30., 60.];

pupil = falco_gen_pupil_ellipse(inputs);

figure(1); imagesc(pupil); axis xy equal tight; colorbar; drawnow;

%% Generate only an ellipse

%--Required Inputs in the structure "inputs"
inputsEllipse.Nbeam = 300.; % max aperture radius in samples
inputsEllipse.Narray = 350; % Number of samples across in square output array
inputsEllipse.radiusX = 0.5; % x-radius of ellipse [pupil diameters]
inputsEllipse.radiusY = 0.2; % y-radius of ellipse [pupil diameters]
inputsEllipse.clockingDegrees = 20; % clocking of the pupil [degrees]

ellipse = falco_gen_ellipse(inputsEllipse);

figure(2); imagesc(ellipse); axis xy equal tight; colorbar; drawnow;
