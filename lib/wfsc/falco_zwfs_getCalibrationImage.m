function Ical = falco_zwfs_getCalibrationImage(mp)
% Simluates the calibration image for the Zernike wavefront sensor using
% FALCO model

    mp.wfs.mask.FPMampFac = 1.0;
    mp.wfs.mask.depth = 0;% Locally modifying the mask depth to zero  
    Ical = falco_zwfs_sim_image(mp);
    
	if(mp.wfs.cam.Npix~=mp.wfs.cam.Narr || any(mp.wfs.cam.centerPixOffset ~=0))
        Ical = fourierPixelate(Ical,mp.wfs.cam.Npix,mp.wfs.cam.centerPixOffset);
    end

end