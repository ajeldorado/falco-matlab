% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear; 

%%--Add to the MATLAB Path
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

input.Npad = 600;
input.Nbeam = 500; % number of points across the pupil diameter
input.OD = 1; % Outer radius (fraction of Nbeam) 
input.ID = 0.2;% Inner radius (zero if you want an off-axis telescope)
input.Nstrut = 4;% Number of struts 
input.angStrut = [0 90 180 270];%Angles of the struts (deg)
input.wStrut = 0.005; % Width of the struts (fraction of pupil diam.)
input.centering = 'pixel';%'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
input.flagPROPER = false;
input.stretch = 1;% Stretch the horizontal axis to create elliptical beam 

PUPIL = falco_gen_pupil_Simple( input );

figure(1);
imagesc(PUPIL);
axis image;
set(gca,'ydir','normal'); colorbar;

%--Check that the centering is correct.
switch input.centering
    case 'interpixel'
        PUPILtemp = PUPIL;
    case 'pixel'
         PUPILtemp = PUPIL(2:end,2:end); %--For pixel centering only
end
figure(2); imagesc(PUPILtemp-rot90(PUPILtemp,2)); axis xy equal tight; colorbar; title('Symmetry Check');
