% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Demonstrate all the options for falco_gen_azimuthal_phase_mask()
%

clear

% Required Inputs
inputs.type = 'staircase';  % Options are 'vortex', 'cos', 'sectors', and 'staircase'.
inputs.N = 200; % pixels across the 
inputs.charge = 4; % charge of the mask (makes most sense for vortex)
inputs.Nsteps = 6; % number of steps per 2*pi radians. For 'staircase' only
inputs.phaseScaleFac = 1; % Factor to apply uniformly to the phase. Used to add chromaticity.

% % Optional Inputs
% inputs.centering = 'pixel';
% inputs.xOffset = 5.5; % [pixels]
% inputs.yOffset = -10; % [pixels]
% inputs.clocking = 0; % [degrees]

% Check staircase
mask = falco_gen_azimuthal_phase_mask(inputs);
figure(1); imagesc(angle(mask)); axis xy equal tight; colorbar; colormap gray; drawnow;

% Check clocking
inputs.clocking = 30; % [degrees]
mask = falco_gen_azimuthal_phase_mask(inputs);
figure(2); imagesc(angle(mask)); axis xy equal tight; colorbar; colormap gray; drawnow;

% Check lateral offsets and rotation together
inputs.clocking = 30; % [degrees]
inputs.xOffset = 30.2; % [pixels]
inputs.yOffset = -14.5; % [pixels]
mask = falco_gen_azimuthal_phase_mask(inputs);
figure(3); imagesc(angle(mask)); axis xy equal tight; colorbar; colormap gray; drawnow;

% Generate staircase with PROPER to get non-binary edges
phase = falco_gen_spiral_staircase(inputs);
mask = exp(1j*phase);
figure(13); imagesc(angle(mask)); axis xy equal tight; colorbar; colormap gray; drawnow;

%%
clear inputs

% Check vortex
inputs.type = 'vortex';  % 'vortex', 'cos', 'sectors', and 'staircase.'
inputs.N = 200; % pixels across the 
inputs.charge = 6; % charge of the mask (makes most sense for vortex)
inputs.phaseScaleFac = 1; % Factor to apply uniformly to the phase. Used to add chromaticity.
mask = falco_gen_azimuthal_phase_mask(inputs);
figure(4); imagesc(angle(mask)); axis xy equal tight; colorbar; colormap gray; drawnow;

% Check sectors
inputs.type = 'sectors';  % 'vortex', 'cos', 'sectors', and 'staircase.'
mask = falco_gen_azimuthal_phase_mask(inputs);
figure(5); imagesc(angle(mask)); axis xy equal tight; colorbar; colormap gray; drawnow;

% Check cos
inputs.type = 'cos';  % 'vortex', 'cos', 'sectors', and 'staircase.'
mask = falco_gen_azimuthal_phase_mask(inputs);
figure(6); imagesc(angle(mask)); axis xy equal tight; colorbar; colormap gray; drawnow;

% Check rotated cos
inputs.clocking = 30; % [degrees]
inputs.type = 'cos';  % 'vortex', 'cos', 'sectors', and 'staircase.'
mask = falco_gen_azimuthal_phase_mask(inputs);
figure(7); imagesc(angle(mask)); axis xy equal tight; colorbar; colormap gray; drawnow;
