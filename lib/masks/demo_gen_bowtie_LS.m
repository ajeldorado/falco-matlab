% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

clear all; 

inputs.Nbeam = 200; % number of points across the incoming beam           
inputs.ID = 0.38; % inner diameter of mask (in pupil diameters)
inputs.OD = 0.91; % outer diameter of mask (in pupil diameters)
inputs.ang = 90; %opening angle of the upper and lower bowtie wedges (degrees)

inputs.clocking = 10;

pupil = falco_gen_bowtie_LS(inputs);


figure(2); imagesc(pupil); axis xy equal tight; colorbar; drawnow;

%%
clear all

clockDeg = 80;%90-60;
centering = 'pixel';

Nbeam = 300;
ID = 0.41;
OD = 0.87;
rocFillet = 0.03;
angDeg = 89;
% clockDeg = 90-60;
upsampleFactor = 11;

DbeamUM = 17.0146e3;
stepSize = 100;
xc = 5e3;
yc = -2234;
fid = [];
aoiDeg = 30;
aoiAxis = 'x';


inputs.Nbeam = Nbeam;     % number of points across the incoming beam           
inputs.ID = ID; % inner diameter of mask (in pupil diameters)
inputs.OD = OD; % outer diameter of mask (in pupil diameters)
inputs.ang = angDeg; %opening angle of the upper and lower bowtie wedges (degrees)
inputs.xShear = xc/DbeamUM;
inputs.yShear = yc/DbeamUM;

inputs.clocking = 30;
maskA = falco_gen_bowtie_LS(inputs);

inputs.Rfillet = 0.03;
maskB = falco_gen_bowtie_LS(inputs);

figure(101); imagesc(maskA); axis xy equal tight; colorbar;
figure(102); imagesc(maskB); axis xy equal tight; colorbar;
figure(103); imagesc(maskA-maskB); axis xy equal tight; colorbar;

