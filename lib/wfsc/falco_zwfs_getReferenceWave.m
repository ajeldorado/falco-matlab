function rw = falco_zwfs_getReferenceWave(mp)
%rw = zwfs_getReferenceWave(mp)
%   Computes the reference wave for the Zernike wavefront sensor using
%   FALCO model

    modvar.wpsbpIndex = 1;
    modvar.sbpIndex = 1;
    
    modvar.whichSource = 'star';
    modvar.lambda = mp.lambda0;
    
%     falcoModelParams_zwfs;
%     mp.FPMampFac = 0.0;
% 	[mp,out] = falco_init_ws([mp.path.config mp.runLabel,'_config.mat']);
%     
    mp.dm1.V = zeros(mp.dm1.Nact);
%     mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad);
    mp.P1.full.E = ones(size(mp.P1.full.E));
    rw = model_ZWFS(mp, modvar, 'refwave');


end

