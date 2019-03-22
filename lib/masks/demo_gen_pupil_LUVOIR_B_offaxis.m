% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear; 
% addpath('segMirrorFunctions');
%%--Add to the MATLAB Path
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

apDia = 500;
input.Nbeam = apDia/0.925; % number of points across the pupil diameter
input.wGap = 6e-3/7.989*apDia; % samples
input.numRings = 4;% Number of rings in hexagonally segmented mirror 
input.Npad = 2^(nextpow2(apDia));
input.ID = 0; % central obscuration radius 
input.OD = 1; % pupil outer diameter, can be < 1
input.Nstrut = 0;% Number of struts 
input.angStrut = [];%Angles of the struts (deg)
input.wStrut = []; % Width of the struts (fraction of pupil diam.)

missingSegments = ones(1,hexSegMirror_numSegments( input.numRings ));
for index = 0:5
    missingSegments(38+index*4) = 0;
end
input.missingSegments = missingSegments;

PUPIL = falco_gen_pupil_customHex( input );

%%

[rows,cols]=size(PUPIL);
[X,Y] = meshgrid(-rows/2:rows/2-1);
[~,RHO] = cart2pol(X,Y);

%%

figure(1);
imagesc(PUPIL+(RHO<apDia/2));
axis image;
set(gca,'ydir','normal');
