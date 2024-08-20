%% dzm_build_q

function Q = dzm_build_q(mp, ev, sbp_texp)

ev.Q = zeros(ev.SS, ev.SS, floor(ev.SL/ev.BS), mp.Nsbp);

for iSubband = 1:mp.Nsbp
    disp(['Assembling Q for subband ', num2str(iSubband)]);
   
    G_reordered = ev.G_tot_drift(:,:,iSubband);
    dm_drift_covariance = eye(size(G_reordered, 2)) * (mp.drift.presumed_dm_std^2);

    for i = 0:floor(ev.SL/ev.BS)-1
        ev.Q(:,:,i+1,iSubband) = G_reordered((i)*ev.BS+1:(i+1)*ev.BS,:) * dm_drift_covariance * G_reordered(i*ev.BS+1:(i+1)*ev.BS,:).' * sbp_texp(iSubband) * (ev.e_scaling(iSubband)^2);
    end
end

Q = ev.Q;

end
