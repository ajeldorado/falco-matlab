function plotParamScanData(mp_arr)

    for i = 1:length(mp_arr)
        mp = mp_arr{i};
        niter = mp.Nitr;
        contrast = zeros(mp_arr{1}.Nitr, 1);
        for ii = 1:niter
            iSubband = ceil(mp.Nsbp / 2);
            I = falco_get_sbp_image(mp, iSubband);
            mean_contrast = mean(I(mp.Fend.score.mask));
            contrast(ii) = mean_contrast;
        end
        figure
        plot(1:niter, contrast);
        out_name = ['Drift_', num2str(mp.drift.magnitude), 'Dither_', num2str(mp.est.dither), '.jpg'];
        saveas(gcf, fullfile(mp.path.ws, out_name))
    end


end