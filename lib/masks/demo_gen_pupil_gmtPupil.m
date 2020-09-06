clear;

%--Library locations

% mp.path.falco  = ;
% mp.path.proper = ;

% addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
% addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


addpath('gmtPupilFunctions');


%% Amplitude only

apDia = 1000;
input.Nbeam = apDia; % number of points across the pupil diameter
input.Npad = 2^(nextpow2(apDia));
input.centralID = 3.2/25.448; % central obscuration fraction
input.segDia = apDia/25.448*8.4; 
input.segID = 0.05; % central obscuration fraction of outer segments 
input.OD = 2; % pupil outer diameter, can be < 1

PUPIL = falco_gen_pupil_gmtPupil( input );

figure(1);
imagesc((PUPIL));
axis image;
set(gca,'ydir','normal');


%% Piston and tip-tilt errors 
piston_std = 1/100; %waves
tilt_std = 1/10; %lambda/D

numSegments = 7;

input.pistons = piston_std.*randn(1,numSegments);% Vector of pistons per segment (waves)
input.tiltxs = tilt_std*randn(1,numSegments); % Vector of x-tilts per segment (waves/apDia)
input.tiltys = tilt_std*randn(1,numSegments);% Vector of y-tilts per segment (waves/apDia)

PUPIL = falco_gen_pupil_gmtPupil( input );

figure(2);
subplot(1,2,1);
imagesc(abs(PUPIL));
axis image;
set(gca,'ydir','normal');
subplot(1,2,2);
imagesc(angle(PUPIL));
axis image;
set(gca,'ydir','normal');


%% Add segment level aberrations 

input.pistons = zeros(1,numSegments);% Vector of pistons per segment (waves)
input.tiltxs = zeros(1,numSegments); % Vector of x-tilts per segment (waves/apDia)
input.tiltys = zeros(1,numSegments);% Vector of y-tilts per segment (waves/apDia)


coeff = 1/10; % segment-level Zernike error in units of waves rms 
noll_indices = 6; % Noll indicies


% Loop over the segments and fill in the loworder_struct 
for segment_index = 1:numSegments
    loworder_struct(segment_index).noll_indices = noll_indices;
    loworder_struct(segment_index).waves_rms = coeff;
end

input.loworder_struct = loworder_struct;

PUPIL = falco_gen_pupil_gmtPupil( input );

figure(3);
subplot(1,2,1);
imagesc(abs(PUPIL));
axis image;
set(gca,'ydir','normal');
subplot(1,2,2);
imagesc(angle(PUPIL));
axis image;
set(gca,'ydir','normal');

%% Add ensemble of segment level aberrations 

coeff_std = 1/50; % standard deviataion of the randomly-generated
                  % segment-level Zernike errors in units of waves rms 
noll_indices = 4:11; % Noll indicies excluding piston, tip, tilt 


% Loop over the segments and fill in the loworder_struct 
for segment_index = 1:numSegments
    loworder_struct(segment_index).noll_indices = noll_indices;
    loworder_struct(segment_index).waves_rms = coeff_std*randn(1,numel(noll_indices));
end

input.loworder_struct = loworder_struct;

PUPIL = falco_gen_pupil_gmtPupil( input );

figure(4);
subplot(1,2,1);
imagesc(abs(PUPIL));
axis image;
set(gca,'ydir','normal');
subplot(1,2,2);
imagesc(angle(PUPIL));
axis image;
set(gca,'ydir','normal');

%% Add ensemble of segment level aberrations with offset 

% Offset the 7th segment by 1% of the pupil diam
input.offset = zeros(numSegments,2);
input.offset(7,1) = 0;
input.offset(7,2) = apDia*0.05;

PUPIL = falco_gen_pupil_gmtPupil( input );

figure(5);
subplot(1,2,1);
imagesc(abs(PUPIL));
axis image;
set(gca,'ydir','normal');
subplot(1,2,2);
imagesc(angle(PUPIL));
axis image;
set(gca,'ydir','normal');

