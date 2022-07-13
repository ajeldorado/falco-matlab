function [mp,ev] = initialize_ekf_maintenance(mp,ev)

ev.peak_psf_counts = zeros(mp.Nsbp);

for iSubband = 1:mp.Nsbp
    if mp.flagSim
        ev.peak_psf_counts(si) = tb.info.sbp_texp(si)*tb.info.PSFpeaks(si);
    else
    
    end

end


end
