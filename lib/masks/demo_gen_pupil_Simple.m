% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear; 

% %%--Add to the MATLAB Path
% mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
% mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
% addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

inputs.Npad = 700;
inputs.Nbeam = 500; % number of points across the pupil diameter
inputs.OD = 1; % Outer radius (fraction of Nbeam) 
inputs.ID = 0.2;% Inner radius (zero if you want an off-axis telescope)
inputs.angStrut = [0 90 180 270];%Angles of the struts (deg)
inputs.wStrut = 0.01; % Width of the struts (fraction of pupil diam.)
inputs.centering = 'pixel';%'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
inputs.stretch = 1;% Stretch the horizontal axis to create elliptical beam 

% Optional Inputs
inputs.clocking = 15; % CCW rotation. Doesn't work with flag HG. [degrees]
inputs.xShear = 0.1;
inputs.yShear = -0.15;
% inputs.flagHG = true; % Use hyper-gaussians for edge antialiasing

pupil = falco_gen_pupil_Simple( inputs );

figure(1);
imagesc(pupil);
axis image;
set(gca,'ydir','normal'); colorbar; colormap gray;

%--Check that the centering is correct.
switch inputs.centering
    case 'interpixel'
        PUPILtemp = pupil;
    case 'pixel'
         PUPILtemp = pupil(2:end,2:end); %--For pixel centering only
end
% figure(2); imagesc(PUPILtemp-rot90(PUPILtemp,2)); axis xy equal tight; colorbar; title('Symmetry Check');


%% Check if it still runs with all optional inputs left out
clear all


inputs.Npad = 600;
inputs.Nbeam = 500; % number of points across the pupil diameter
inputs.OD = 1; % Outer radius (fraction of Nbeam) 
% inputs.ID = 0.2;% Inner radius (zero if you want an off-axis telescope)
% inputs.Nstrut = 4;% Number of struts 
% inputs.angStrut = [0 90 180 270];%Angles of the struts (deg)
% inputs.wStrut = 0.005; % Width of the struts (fraction of pupil diam.)
inputs.centering = 'pixel';%'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
% inputs.flagPROPER = false;
% inputs.stretch = 1;% Stretch the horizontal axis to create elliptical beam 

pupil = falco_gen_pupil_Simple(inputs);

figure(3); imagesc(pupil); axis xy equal tight; colorbar; colormap gray; title('Simplest Pupil');
