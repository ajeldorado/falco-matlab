% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

clear all; 

inputs.Nbeam = 200; % number of points across the incoming beam           
inputs.ID = 0.38; % inner diameter of mask (in pupil diameters)
inputs.OD = 0.91; % outer diameter of mask (in pupil diameters)
inputs.ang = 90; %opening angle of the upper and lower bowtie wedges (degrees)

inputs.clocking = 10;

pupil = falco_gen_bowtie_LS(inputs);


figure(2); imagesc(pupil); axis xy equal tight; colorbar; drawnow;

% %% Generate only an ellipse
% 
% %--Required Inputs in the structure "inputs"
% inputsEllipse.Nbeam = 300.; % max aperture radius in samples
% inputsEllipse.Narray = 350; % Number of samples across in square output array
% inputsEllipse.radiusX = 0.5; % x-radius of ellipse [pupil diameters]
% inputsEllipse.radiusY = 0.2; % y-radius of ellipse [pupil diameters]
% inputsEllipse.clockingDegrees = 20; % clocking of the pupil [degrees]
% 
% ellipse = falco_gen_ellipse(inputsEllipse);
% 
% figure(2); imagesc(ellipse); axis xy equal tight; colorbar; drawnow;
