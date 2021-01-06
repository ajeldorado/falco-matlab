% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the final LUVOIR A telescope pupil in Matlab using PROPER
%
%--Coordinates of hex segments to skip:
% 1 13 114 115 126 127
% 1 12 113 114 125 126

clear all

% %%--Add to the MATLAB Path
% mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
% mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
% addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

inputs.Nbeam = 1000;
inputs.magFac = 1;
inputs.wStrut = 1/100;
% inputs.strutWidth_m = 0; % Optional input, strutwidth in meters 

pupil = falco_gen_pupil_LUVOIR_A_final(inputs);

figure(1); imagesc(pupil); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar; drawnow;

figure(2); imagesc(pupil(2:end,2:end)-fliplr(pupil(2:end,2:end))); axis xy equal tight; title('Symmetry Check: Differencing','Fontsize',20); colorbar; drawnow;

inputs.magFac = 0.3;
inputs.clock_deg = 30;
inputs.xShear = 0.4;
pupil2 = falco_gen_pupil_LUVOIR_A_final(inputs);

figure(2); imagesc(pupil2); axis xy equal tight; colorbar; drawnow;
