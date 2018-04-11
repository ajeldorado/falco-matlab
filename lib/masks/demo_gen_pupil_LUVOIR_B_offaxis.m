% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear; 
addpath('segMirrorFunctions');


apDia = 500;
input.Nbeam = apDia/0.925; % number of points across the pupil diameter
input.gapWidth = 6e-3/7.989*apDia; % samples
input.numRings = 4;% Number of rings in hexagonally segmented mirror 
input.Npad = 2^(nextpow2(apDia));
input.ID = 0; % central obscuration radius 
input.OD = 1; % pupil outer diameter, can be < 1
input.num_strut = 0;% Number of struts 
input.strut_angs = [];%Angles of the struts (deg)
input.strut_width = []; % Width of the struts (fraction of pupil diam.)

missingSegments = ones(1,hexSegMirror_numSegments( input.numRings ));
for index = 0:5
    missingSegments(38+index*4) = 0;
end
input.missingSegments = missingSegments;

PUPIL = gen_pupil_customHex( input );

%%

[rows,cols]=size(PUPIL);
[X,Y] = meshgrid(-rows/2:rows/2-1);
[~,RHO] = cart2pol(X,Y);

%%

figure(1);
imagesc(PUPIL+(RHO<apDia/2));
axis image;
set(gca,'ydir','normal');
