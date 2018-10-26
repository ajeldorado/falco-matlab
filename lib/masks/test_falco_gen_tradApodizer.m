clear; close all;

%% Make a pupil 
apRad = 500; % radius in samples 
N = 2^nextpow2(2*apRad);% size of array 
Npad = 4*N; %size of padded array 

[X,Y] = meshgrid(-N/2:N/2-1);
[~,RHO] = cart2pol(X,Y);
xvals = X(1,:); yvals = Y(:,1);

EP = exp(-(RHO/apRad).^5000);
PSF0 = abs(fftshift(fft2(fftshift(padOrCropEven(EP,Npad))))).^2;
normI = max(max(PSF0));
PSF0 = PSF0/normI; 

figure;
imagesc(xvals,yvals,EP);
axis image; 
colorbar; 

%% Compute the apodizer
APOD = falco_gen_tradApodizer(EP,2*apRad,3,10,false);

figure;
imagesc(xvals,yvals,APOD);
axis image; 
colorbar; 

%% Compute the PSF

APOD = padOrCropEven(APOD,4*N);

FP = fftshift(fft2(fftshift(APOD)));
PSF = abs(FP).^2/normI;
lamOverD = Npad / (2*apRad);

xvals_fp = -Npad/2:Npad/2-1;
yvals_fp = xvals_fp';

figure;
imagesc(xvals_fp/lamOverD,yvals_fp/lamOverD,log10(PSF0));
axis image; 
colorbar; caxis([-8 0])
axis([-10 10 -10 10]);

figure;
imagesc(xvals_fp/lamOverD,yvals_fp/lamOverD,log10(PSF));
axis image; 
colorbar; caxis([-8 0])
axis([-10 10 -10 10]);