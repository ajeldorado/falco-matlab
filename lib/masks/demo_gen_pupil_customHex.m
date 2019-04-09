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
% addpath('segMirrorFunctions');

%% Amplitude only (Keck pupil)

apDia = 1000;
input.Nbeam = apDia; % number of points across the pupil diameter
input.wGap = 25.4/10918*apDia/2; % samples
input.numRings = 3;% Number of rings in hexagonally segmented mirror 
input.Npad = 2^(nextpow2(apDia));
input.ID = 2600/10949; % central obscuration radius 
input.OD = 2; % pupil outer diameter, can be < 1
input.Nstrut = 6;% Number of struts 
input.angStrut = [0 60 120 180 240 300];%Angles of the struts (deg)
input.wStrut = 25.4/10949; % Width of the struts (fraction of pupil diam.)

PUPIL = falco_gen_pupil_customHex( input );

figure(1);
imagesc((PUPIL));
axis image;
set(gca,'ydir','normal');

%% Piston and tip-tilt errors (Keck pupil)
piston_std = 1/100; %waves
tilt_std = 1/10; %lambda/D

numSegments = hexSegMirror_numSegments( input.numRings );

input.pistons = piston_std.*randn(1,numSegments);% Vector of pistons per segment (waves)
input.tiltxs = tilt_std*randn(1,numSegments); % Vector of x-tilts per segment (waves/apDia)
input.tiltys = tilt_std*randn(1,numSegments);% Vector of y-tilts per segment (waves/apDia)

PUPIL = falco_gen_pupil_customHex( input );

figure(2);
subplot(1,2,1);
imagesc(abs(PUPIL));
axis image;
set(gca,'ydir','normal');
subplot(1,2,2);
imagesc(angle(PUPIL));
axis image;
set(gca,'ydir','normal');

%% LUVOIR A

input.Nbeam = apDia; % number of points across the pupil diameter
input.wGap = 25.4/10918*apDia/2; % samples
input.numRings = 7;% Number of rings in hexagonally segmented mirror 
input.Npad = 2^(nextpow2(apDia));
input.ID = 0; % central obscuration radius 
input.OD = 2; % pupil outer diameter, can be < 1
input.Nstrut = 0;% Number of struts 
input.angStrut = [];%Angles of the struts (deg)
input.wStrut = 0; % Width of the struts (fraction of pupil diam.)

numSegments = hexSegMirror_numSegments( input.numRings );

input.pistons = piston_std.*randn(1,numSegments);% Vector of pistons per segment (waves)
input.tiltxs = tilt_std*randn(1,numSegments); % Vector of x-tilts per segment (waves/apDia)
input.tiltys = tilt_std*randn(1,numSegments);% Vector of y-tilts per segment (waves/apDia)

numSegments = hexSegMirror_numSegments( input.numRings );

PUPIL = falco_gen_pupil_customHex( input );

figure(3);
subplot(1,2,1);
imagesc(abs(PUPIL));
axis image;
set(gca,'ydir','normal');
subplot(1,2,2);
imagesc(angle(PUPIL));
axis image;
set(gca,'ydir','normal');
