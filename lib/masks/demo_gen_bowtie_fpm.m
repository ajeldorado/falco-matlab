% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

clear all; 

inputs.pixresFPM = 6; %--pixels per lambda_c/D
inputs.rhoInner = 2.6; % radius of inner FPM amplitude spot (in lambda_c/D)
inputs.rhoOuter = 9.4; % radius of outer opaque FPM ring (in lambda_c/D)
inputs.ang = 60 ;
inputs.centering = 'pixel';

%--Optional Inputs
inputs.xOffset = 5.5;
inputs.yOffset = -10;
% inputs.Narray = 512;

fpm = falco_gen_bowtie_FPM(inputs);


figure(1); imagesc(fpm); axis xy equal tight; colorbar; drawnow;
