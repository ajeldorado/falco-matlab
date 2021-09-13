function Isum = falco_zwfs_sim_image(mp)

    %--Compute the DM surfaces outside the full model to save lots of time
    if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
    if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end

    Isum = 0; % Initialize image

    for si=1:mp.wfs.Nsbp
        modvar.sbpIndex = si;

        for wi=1:mp.wfs.Nwpsbp
            modvar.whichSource = 'star';
            modvar.wpsbpIndex = wi;
            Etemp = model_ZWFS(mp, modvar);

            Isum = Isum + (abs(Etemp).^2)*mp.wfs.sbp_weights(si)*mp.wfs.lambda_weights(wi);
        end 

    end
    
	if(mp.wfs.cam.Npix~=mp.wfs.cam.Narr || any(mp.wfs.cam.centerPixOffset ~=0))
        Isum = fourierPixelate(Isum,mp.wfs.cam.Npix,mp.wfs.cam.centerPixOffset);
    end

end
