% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear; 
addpath('../utils');

Nbeam = 500;
Npad = 512;
PUPIL = falco_gen_pupil_dst_LUVOIR_B( Nbeam, Npad );

%%

[rows,cols]=size(PUPIL);
[X,Y] = meshgrid(-rows/2:rows/2-1);
[~,RHO] = cart2pol(X,Y);

%%

figure(1);
imagesc(PUPIL+(RHO<Nbeam/2));
axis image;
set(gca,'ydir','normal');
