%% Read in Lorenzo Metasurface Maps
%Run this script before simplefpmtester to load lorenzo's votcube bc it's 
%too big to save and read in!

close all
clear all

%design = "vortex"; %
design="roddier";
trans  = "var"; %"const";
lams   = [3.4,3.6,3.8,4.,4.2];
% i_lam  = 2;%2;%3;%4;%5
vortcube = [];%zeros(4);


% mask = exp(phaseScaleFac*1j*roddiervort);


for i_lam = 1:length(lams)

    Npad   = 4096; %4096; %8192; % %2^13;
    charge = 6;
    theta  = angle(falco_gen_vortex_mask(1,Npad)); % gen_v_mask  with charge 1 means having an array with the azi angle for each pxl
    z1 = 2.40483; % wiki Bessel zeros
    z2 = 5.52008; % using z2 or z3 results in more than 2pi phase coverage needed
    z3 = 8.65373; %
    
    %create roddier+sawtooth phase ideal pattern
    coords = generateCoordinates(Npad);% Creates NxN arrays with coordinates 
    roddiervort = 0.* coords.THETA;
    domain = (coords.THETA >= 0);
    roddiervort(domain) = charge*rem(coords.THETA(domain),2*pi./charge);
    domain = (coords.THETA >= -pi) & (coords.THETA < 0);
    roddiervort(domain) = charge*rem((coords.THETA(domain)+pi),2*pi./charge);
    res = 1024;
    R1 = (coords.RHO <= 0.53*res);
    roddiervort(R1) =roddiervort(R1) + 0.5*2*pi;
    
    %TRANSLATE TO METASURFACE MAP
    if design=="vortex", [pss,Tss,sss] = b30_library_siz_auto(angle(exp(1i.*charge.*theta)),i_lam); end % VORTEX: add wavelength dependent phase. 
    if design=="roddier", [pss,Tss,sss] = b30_library_siz_auto(angle(exp(1i.*roddiervort)),i_lam); end % RODDIER: add wavelength dependent phase. 
    if design=="cosine", [pss,Tss,sss] = b30_library_siz_auto(angle(exp(1i.*z1.*cos(charge.*theta)))+z1-pi,i_lam); end % COSINE: add wavelength dependent phase. 
    vortex1 = exp(1i.*pss); % phase only
    if trans=="const", vortex = vortex1; end
    if trans=="var", vortex = vortex1.*sqrt(Tss); end % add phase dependent transmission values
    
    vortcube = cat(3,vortcube,vortex); %(:,:,i_lam);
    
end


%% 1D profile

%%% for (sawtooth) MS vortex
lp=6; % topological charge
N=200 % azimuthal sampling
[phase,transmission,size]=b30_library_siz_auto(angle(exp(1i.*linspace(0,lp*2*pi,N))),linspace(1,5,5));
i_lam=5; % wavelength index 3 is "perfect" (design wavelength)
offsetphase = phase(:,i_lam);
figure()
plot(offsetphase,'o-')
figure()
plot(phase(:,i_lam)-phase(:,3),'o-')

%%Fourier Transform:
fftmetasurface = abs(fftshift(fft(offsetphase)))/N;

%%Frequency specifications:
Fs = N;
dF = Fs/N;  %1                  % hertz
f = -Fs/2:dF:Fs/2-dF;           % hertz

%%Plot the spectrum:
figure(311);
hold on
% plot(f,fftmetasurface);
stem(f,fftmetasurface,'d','MarkerSize',10,'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]);
xlim([-2 20]);
ylim([1E-5 10])
xticks([0:4:18])
set(gca, 'YScale', 'log')
ylabel('|C_m|^2'); 
xlabel('Mode');
title("FFT of Metasurface Phase Profile")



%%
save lorenzoroddierattempt.mat vortcube lams design
% save lorenzovortexmap3.mat vortex



%% Read in Skyler Metasurface Maps


SkylerMapsDir = '/Users/niyatid/falco-matlab/metasurfaces/';
% Input parameters 

%-- Define Filenames
  % NOTE::: Notice the "%d" replacing what should be the input wavelengths
  % I assume all files will have the same format so this %d lets us
  % programatically load the files while also knowing the input wavelengths
  % Input matrices should be square, powers of 2 (256x256, 2048x2048, 4096x4096, etc.)
  % TLDR ==>> only provide one filename for each type of file and use %d
  %           to replace wavelength values and only provide square matrices
ampflnm = fullfile(SkylerMapsDir,'Transmission_%2dum_SiSi3N4.txt');     % Amplitude (transmission) file
phsflnm = fullfile(SkylerMapsDir,'phase_%2dum_SiSi3N4.txt');     % Amplitude (transmission) file

%-- Define spatial sampling in pupil plane
  % (This should be your input grid size: 256x256, 1024x1024, 4096x4096, etc.)
  % NOTE::: powers of 2 make the code run significantly faster due to optimizations in the backend
N = 2^8;

%-- Define wavelength simulation info
% Provide list of avalable input file wavelengths
  % These should only include the wavelengths for which there are files
  % available to load
lams = [544, 598, 645, 700];     % [microns]
% Now define number of wavelenth samples to use for the sim:
  % if != to length of inlambdas, will linearly interpolate vortex data
  % NOTE::: strongly recommend ODD NUMBER and at least as many as inlambdas
  % Look at output print statements to make sure sampling makes sense...
numWavelengths = 4; % number of discrete wavelengths  to use

%% Load the vortex specification files (If necessary)

% Preallocate matrices to hold vortex data
EPMP = nan(N, N, length(lams));      % vortex phase matrix
EPMA = nan(N, N, length(lams));      % vortex amplitude matrix

% Iterate through provided input wavelengths and load files for each one
for ch = 1:length(lams)
    lam = lams(ch);    % Current wavelength 
    % Load transmission file as amplitude directly
    EPMA(:,:,ch) = readmatrix(sprintf(ampflnm,round(lam)));
    % Load phase file as real-valued element that will go in exponential
    EPMP(:,:,ch) = readmatrix(sprintf(phsflnm,round(lam)));

    % Correct amplitude to be at most 1 (some files showed >1 transmission which is unphysical)
    EPMA = min(EPMA, ones(size(EPMA)));

end


EPMP = exp(1i*EPMP);
% Combine vortex phase and amplitude into single mask
vortcube = EPMA.*EPMP;




function coords = generateCoordinates( N )
    %[ X,Y,THETA,RHO,xvals,yvals ] = generateCoordinates( N )
    %   Generates sample centered coordinate system (both cartesian and polar)

        % Create coordinate system 
        [X,Y] = meshgrid(-N/2:N/2-1);
        [THETA,RHO] = cart2pol(X,Y);
        xvals = X(1,:);yvals = Y(:,1);

        coords.N = N;
        coords.X = X;
        coords.Y = Y;
        coords.THETA = THETA;
        coords.RHO = RHO;
        coords.xvals = xvals;
        coords.yvals = yvals;
end