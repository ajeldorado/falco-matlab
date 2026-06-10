function [mp, ev] = falco_take_mono_images(mp, ev)

    Nwpsbp_old = mp.Nwpsbp;
    lambda0_old = mp.lambda0;
    sbp_texp_old = mp.tb.info.sbp_texp;

    wl_bound_low = mp.lambda0 - mp.lambda0 * fracBW_old/2;
    wl_bound_high = mp.lambda0 + mp.lambda0 * fracBW_old/2;
    fracBW_new = 0.025;
    wl_cent_low = wl_bound_low/(1-fracBW_new/2);
    wl_cent_high = wl_bound_high/(1+fracBW_new/2);
    mp.fracBW = fracBW_new;
    mp.Nwpsbp = 1;
    mp.tb.info.sbp_texp = 60;
    images_mono = zeros(500, 500, 3, mp.Nsbp);
    bands = [wl_cent_low, lambda0_old, wl_cent_high];
    for iSubband = 1:mp.Nsbp
        for b = 1:length(bands)
            mp.lambda0 = bands(b);
            mp.sbp_centers = mp.lambda0;
            images_mono(:,:,b, iSubband) = falco_get_sbp_image(mp, iSubband);
            if ~mp.flagSim
                images_mono(:,:,b, iSubband) = images_mono(:,:,b, iSubband) * (mp.tb.info.PSFpeaks(1) * mp.tb.info.sbp_texp) / (tb.info.PSFpeaks_mono(b) * mp.tb.info.sbp_texp);
            end
        end
    end
    fitswrite(images_mono, fullfile(mp.path.ws,sprintf("monochrome_images_itr%d.fits", ev.Itr)));
    mp.fracBW = fracBW_old;
    mp.Nwpsbp = Nwpsbp_old;
    mp.tb.info.sbp_texp = sbp_texp_old;
    mp.lambda0 = lambda0_old;
    mp.sbp_centers = mp.lambda0;
end