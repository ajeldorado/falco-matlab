clear ;

I_in = checkerboard(100);
Npix = 64;
centerPixOffset = [0 0];%(row,col) 
% centerPixOffset = [-0.5 -0.5];%(row,col) 

%%

%Icam = fourierPixelate(I_in,Npix);
Icam = fourierPixelate(I_in,Npix,centerPixOffset);

figure(1);
imagesc(I_in);
axis image;
colorbar;

figure(2);
imagesc(Icam);
axis image;
colorbar;
