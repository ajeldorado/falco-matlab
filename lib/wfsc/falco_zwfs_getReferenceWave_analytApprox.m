function rw = falco_zwfs_getReferenceWave_analytApprox(mp)
% Computes the reference wave for the Zernike wavefront sensor using
% an analytical approximation (2nd order Taylor expansion)

    dhat = 2*mp.wfs.mask.Rin;   % FPM diameter in units of lam0/D
    Nbeam = mp.wfs.cam.Nbeam;   % Number of pixels across the beam at WFS camera
    Narr = mp.wfs.cam.Narr;     % Number of pixels across the full WFS camera image
    
    [X,Y] = meshgrid(-Narr/2:Narr/2-1);
    RHO = sqrt(X.^2+Y.^2);
    rw = 1 - besselj(0,pi*dhat/2) ...
        - (pi*dhat/2)^2*besselj(2,pi*dhat/2).*(RHO/Nbeam).^2;

end