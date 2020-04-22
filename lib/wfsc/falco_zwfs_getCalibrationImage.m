function Ical = falco_zwfs_getCalibrationImage(mp)
% Simluates the calibration image for the Zernike wavefront sensor using
% FALCO model

    mp.wfs.mask.FPMampFac = 1.0;
    mp.wfs.mask.depth = 0;% Locally modifying the mask depth to zero  
    Ical = falco_zwfs_sim_image(mp);

end