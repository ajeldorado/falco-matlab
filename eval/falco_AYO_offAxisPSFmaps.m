%--Throughput part:
SepVec = (0:mp.Fend.Neta/2-1)/mp.Fend.res; %--Vector of angular separations at which to compute the core throughput (lambda_central/D). Points outside the dark hole are removed automatically.

offaxisPSFcube = zeros(mp.Fend.Neta,mp.Fend.Nxi,length(SepVec));

thput_curve = zeros(length(SepVec),1);
transList = zeros(length(SepVec),1);
for ti = 1:length(SepVec)

    mp.thput_eval_x = SepVec(ti);
    mp.thput_eval_y = 0;

    [XIS,ETAS] = meshgrid(mp.Fend.xisDL - mp.thput_eval_x, mp.Fend.etasDL - mp.thput_eval_y);
    mp.Fend.RHOS = sqrt(XIS.^2 + ETAS.^2);
    mp.maskHMcore = 0*mp.Fend.RHOS;
    mp.maskCore  = 0*mp.Fend.RHOS;
    mp.maskCore(mp.Fend.RHOS<=mp.thput_radius) = 1;
    
    modvar.x_offset = mp.thput_eval_x;
    modvar.y_offset = mp.thput_eval_y;
    modvar.sbpIndex = mp.si_ref; 
    modvar.wpsbpIndex = mp.wi_ref;
    modvar.whichSource = 'offaxis';
    ImTemp = falco_sim_image(mp, modvar);
%     if(mp.flagPlot); figure(324); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,ImTemp/mp.sumPupil*mp.Fend.full.I00(mp.si_ref,mp.wi_ref)); axis xy equal tight; title('Off-axis PSF for Throughput Calculation','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end

    offaxisPSFcube(:,:,ti) = ImTemp/mp.sumPupil*mp.Fend.full.I00(mp.si_ref,mp.wi_ref);
    transList(ti) = sum(sum(offaxisPSFcube(:,:,ti)));
    
% maskHM = 0*mp.Fend.RHOS;
% maskH264M(ImTemp>=1/2*max(max(ImTemp))) = 1;
% mp.maskHMcore = maskHM.*mp.maskCore;
mp.maskHMcore = mp.maskCore; %--Do not block out points below half max

%figure(325); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,mp.maskCore); axis xy equal tight; drawnow;
thput_curve(ti) = sum(ImTemp(mp.maskHMcore==1))/mp.sumPupil*mp.Fend.full.I00(mp.si_ref,mp.wi_ref);

end

if(mp.flagPlot)
    figure(10); plot(SepVec,thput_curve,'-bo','Linewidth',2','MarkerSize',6);
    xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
    ylabel('Core Throughput','FontSize',16,'Interpreter','LaTeX'); 
    set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
end

pp = spline(SepVec*mp.Fend.res,transList);

[Xfp,Yfp] = meshgrid((-mp.Fend.Neta/2):(mp.Fend.Neta/2-1));
[~,RHOfp] = cart2pol(Xfp,Yfp);

transMap_fit = ppval(pp,RHOfp);
transMap_fit = transMap_fit.*(RHOfp<=mp.Fend.Neta/2-1);

%% Output to files 
offaxisPSFcube_flnm = [Dir,'offax_psf.fits'];
flnm4 = [Dir,'offax_psf_offset_list.fits'];
transMap_flnm = [Dir,'sky_trans.fits'];

savepath
restoredefaultpath
rehash toolboxcache

import matlab.io.*
try
    fptr = fits.createFile(offaxisPSFcube_flnm);
catch 
    delete(offaxisPSFcube_flnm);
    fptr = fits.createFile(offaxisPSFcube_flnm);
end
fits.createImg(fptr,'double_img',[mp.Fend.Neta,mp.Fend.Nxi,numel(SepVec)]);
fits.writeImg(fptr,offaxisPSFcube);
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
fitsdisp(offaxisPSFcube_flnm,'mode','full');


psfOffsetList = [SepVec',zeros(numel(SepVec),1)];
% psfOffsetList = [SepVec',(0.5/mp.Fend.res)*ones(numel(SepVec),1)];
try
    fptr = fits.createFile(flnm4);
catch 
    delete(flnm4);
    fptr = fits.createFile(flnm4);
end
fits.createImg(fptr,'double_img',[numel(SepVec) 2]);
fits.writeImg(fptr,psfOffsetList);
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
fitsdisp(flnm4,'mode','full');


try
    fptr = fits.createFile(transMap_flnm);
catch 
    delete(transMap_flnm);
    fptr = fits.createFile(transMap_flnm);
end
fits.createImg(fptr,'double_img',[mp.Fend.Neta,mp.Fend.Nxi]);
fits.writeImg(fptr,transMap_fit);
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
fitsdisp(transMap_flnm,'mode','full');

path(pathdef);
