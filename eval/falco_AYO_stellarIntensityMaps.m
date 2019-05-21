% %--NEED TO HAVE MP AS 2ND OUTPUT OF falco_wfsc_loop TO HAVE MP TO USE
% [out,mp] = falco_wfsc_loop(mp);

% mp.full.TipTiltNacross = 5; %numOfPts = 5; %--Radial number of points across to use in the tip/tilt offset grid.
% dia_ld_list = [0 0.005 0.5]; % List of stellar sizes in lam/D
mp.full.TipTiltNacross = 9; %numOfPts = 9; %--Radial number of points across to use in the tip/tilt offset grid.
dia_ld_list = [0:0.005:0.1,0.15:0.05:0.5]; % List of stellar sizes  [lam/D]

mp.lambda0 = 500e-9; %--central wavelength [meters]
mp.full.Dtel = 15; %--diameter of the primary mirror [meters]
mp.full.TTrms = 1; %--RMS pointing jitter [mas]

lambdaOverD_to_mas = mp.lambda0/mp.full.Dtel*180/pi*3600*1000;
Dstar_list = dia_ld_list*lambdaOverD_to_mas;

stellarPSFcube = zeros(mp.Fend.Neta,mp.Fend.Nxi,numel(dia_ld_list));
count = 1;
for Dstar = Dstar_list
    mp.full.Dstar = Dstar;
    fprintf('\nComputing the stellar psf for Dstar = %s lam/D',num2str(Dstar));
    
    Imean = falco_get_summed_image_TipTiltPol(mp);
    
    stellarPSFcube(:,:,count) = Imean/mp.sumPupil*mean(mp.Fend.full.I00(:));

    count = count + 1; % ah ah ah
end

%%
% disp('Shifting the stellar PSF by half a pixel to match AYO standard');
% for index = 1:numel(Dstar_list)
%    
%     stellarPSFcube(:,:,index) = imtranslate(stellarPSFcube(:,:,index),[-0.5 -0.5]);
%     
% end

%%

iPSFstarBB_flnm = [Dir,'stellar_intens.fits'];
flnm2 = [Dir,'stellar_intens_diam_list.fits'];

savepath
restoredefaultpath
rehash toolboxcache

import matlab.io.*
try
    fptr = fits.createFile(iPSFstarBB_flnm);
catch 
    delete(iPSFstarBB_flnm);
    fptr = fits.createFile(iPSFstarBB_flnm);
end
fits.createImg(fptr,'double_img',[mp.Fend.Neta,mp.Fend.Nxi,numel(dia_ld_list)]);
fits.writeImg(fptr,stellarPSFcube);
fits.writeKey(fptr,'PIXSCALE',1/mp.Fend.res,'sample spacing in units of lambda_0/D');
fits.writeKey(fptr,'LAMBDA',lambda0*1e6,'lambda_0, central wavelength (um)');
fits.writeKey(fptr,'MINLAM',mp.lam_array(1)*1e6,'shortest wavelength (um)');
fits.writeKey(fptr,'MAXLAM',mp.lam_array(end)*1e6,'longest wavelength (um)');
fits.writeKey(fptr,'D',Dtel,'Aperture Diameter (m)');
fits.writeKey(fptr,'OBSCURED',obscuredFraction,'Fraction of aperture area obscured');
fits.writeKey(fptr,'XCENTER',mp.Fend.Neta/2+0.5,'the coordinate of the image center');
fits.writeKey(fptr,'YCENTER',mp.Fend.Nxi/2+0.5,'the coordinate of the image center');
fits.writeKey(fptr,'JITTER',TTrms,'RMS jitter per axis in mas');
fits.writeKey(fptr,'UNITS','unitless','units of the values in the array');
fits.closeFile(fptr);
fitsdisp(iPSFstarBB_flnm,'mode','full');

try
    fptr = fits.createFile(flnm2);
catch 
    delete(flnm2);
    fptr = fits.createFile(flnm2);
end
fits.createImg(fptr,'double_img',[numel(dia_ld_list) 1]);
fits.writeImg(fptr,dia_ld_list);
fits.writeKey(fptr,'PIXSCALE',1/mp.Fend.res,'sample spacing in units of lambda_0/D');
fits.writeKey(fptr,'LAMBDA',lambda0*1e6,'lambda_0, central wavelength (um)');
fits.writeKey(fptr,'MINLAM',mp.lam_array(1)*1e6,'shortest wavelength (um)');
fits.writeKey(fptr,'MAXLAM',mp.lam_array(end)*1e6,'longest wavelength (um)');
fits.writeKey(fptr,'D',Dtel,'Aperture Diameter (m)');
fits.writeKey(fptr,'OBSCURED',obscuredFraction,'Fraction of aperture area obscured');
fits.writeKey(fptr,'XCENTER',mp.Fend.Neta/2+0.5,'the coordinate of the image center');
fits.writeKey(fptr,'YCENTER',mp.Fend.Nxi/2+0.5,'the coordinate of the image center');
fits.writeKey(fptr,'JITTER',TTrms,'RMS jitter per axis in mas');
fits.writeKey(fptr,'UNITS','LAMBDA/D','units of the values in the array');
fits.closeFile(fptr);
fitsdisp(flnm2,'mode','full');

path(pathdef);