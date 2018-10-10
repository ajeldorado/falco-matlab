% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear; 

input.Npad = 600;
input.Nbeam = 500; % number of points across the pupil diameter
input.OD = 1; % Outer radius (fraction of Nbeam) 
input.ID = 0.2;% Inner radius (zero if you want an off-axis telescope)
input.num_strut = 4;% Number of struts 
input.strut_angs = [0 90 180 270];%Angles of the struts (deg)
input.strut_width = 0.005; % Width of the struts (fraction of pupil diam.)
input.centering = 'interpixel'; %'interpixel' or 'pixel'; 'pixel' is default
input.flagPROPER = true;

PUPIL = falco_gen_pupil_Simple( input );

figure(1);
imagesc(PUPIL);
axis image;
set(gca,'ydir','normal');

%--Check that the centering is correct.
PUPILtemp = PUPIL(2:end,2:end); %--For pixel centering only
figure(2); imagesc(PUPILtemp-rot90(PUPILtemp,2)); axis xy equal tight; colorbar;