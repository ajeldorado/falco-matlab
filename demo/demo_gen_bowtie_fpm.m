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

%% Compare fillet and no-fillet versions
clear all

inputs.clocking = 0;

inputs.pixresFPM = 6; %--pixels per lambda_c/D
inputs.rhoInner = 2.6; % radius of inner FPM amplitude spot (in lambda_c/D)
inputs.rhoOuter = 9.4; % radius of outer opaque FPM ring (in lambda_c/D)
inputs.ang = 65 ;
inputs.centering = 'pixel';


inputs.xOffset = 2;
inputs.yOffset = -5;

maskA = falco_gen_bowtie_FPM(inputs);

inputs.Rfillet = 0.50;
maskB = falco_gen_bowtie_FPM(inputs);

figure(101); imagesc(maskA); axis xy equal tight; colorbar;
figure(102); imagesc(maskB); axis xy equal tight; colorbar;
figure(103); imagesc(maskA-maskB); axis xy equal tight; colorbar;



%%
% clear all
% 
% clockDeg = 0;
% 
% rhoInner = 2.6;
% rhoOuter = 9.4;
% rocFillet = 0.5;
% pixresFPM = 10;
% angDeg = 65;
% upsampleFactor = 101;
% centering = 'pixel';
% 
% maskA = falco_gen_rounded_bowtie_FPM(rhoInner, rhoOuter, rocFillet, pixresFPM, angDeg, clockDeg, upsampleFactor, centering);
% 
% maskB = falco_gen_rounded_bowtie_FPM_no_rotate(rhoInner, rhoOuter, rocFillet, pixresFPM, angDeg, clockDeg, upsampleFactor, centering);
% 
% figure(101); imagesc(maskA); axis xy equal tight; colorbar;
% figure(102); imagesc(maskB); axis xy equal tight; colorbar;
% figure(103); imagesc(maskA-maskB); axis xy equal tight; colorbar;
