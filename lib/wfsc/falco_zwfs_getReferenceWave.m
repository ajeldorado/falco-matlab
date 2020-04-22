function rw = falco_zwfs_getReferenceWave(mp)
% Computes the reference wave for the Zernike wavefront sensor using
% FALCO model

    modvar.wpsbpIndex = 1;
    modvar.sbpIndex = 1;
    
    modvar.whichSource = 'star';
    modvar.lambda = mp.wfs.lambda0;
   
    % assuming the the DMs and nominal phase are flat
    mp.dm1.V = zeros(mp.dm1.Nact);
    mp.dm2.V = zeros(mp.dm2.Nact);
    mp.P1.full.E = ones(size(mp.P1.full.E));
    rw = model_ZWFS(mp, modvar, 'refwave');
    rw = abs(rw); 

end