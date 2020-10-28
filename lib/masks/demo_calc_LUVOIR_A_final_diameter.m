% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to determine the outer diameter of the LUVOIR-A pupil.

clear all

%%--Add to the MATLAB Path
% mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
% mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
% addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

fac = 2;
Npad = fac*1050;
Nfinal = ceil_even(fac*1049);

Nbeam0 = fac*1000;
inputs.Nbeam = Nbeam0;
inputs.magFac = 1;
inputs.clock_deg = 0;

pupil = falco_gen_pupil_LUVOIR_A_final(inputs);
pupil = pad_crop(pupil, Nfinal);

figure(2); imagesc(pupil); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;

inputs.Npad = Npad; %Number of samples in NxN grid 
inputs.OD = 1; % pupil outer diameter, can be < 1
inputs.Nbeam = Nbeam0*1.000; %1.01824 ;
circle = pad_crop(falco_gen_pupil_Simple(inputs), Nfinal);

figure(1); imagesc(pupil); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;
figure(2); imagesc(circle - pupil); axis xy equal tight; colorbar; drawnow;
fprintf('Min value %.4f for ratio = %.5f\n', min(min(circle - pupil)),inputs.Nbeam/Nbeam0)


%% See if changing segment gaps changes OD --> Not much because dx changes within the function.
clear all;

inputs.Nbeam = 1000;
Nfinal = 1020;

pupil1 = falco_gen_pupil_LUVOIR_A_final(inputs);
pupil1 = pad_crop(pupil1, Nfinal);

figure(1); imagesc(pupil1); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;


inputs.wGap_m = 0;
pupil2 = falco_gen_pupil_LUVOIR_A_final(inputs);
pupil2 = pad_crop(pupil2, Nfinal);

figure(2); imagesc(pupil2-pupil1); axis xy equal tight; title('Difference','Fontsize',20); colorbar;
