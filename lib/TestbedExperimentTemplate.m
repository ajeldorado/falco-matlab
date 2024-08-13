% Exposure times
function tb = TestbedExperimentTemplate(mp)
    tb.info.sbp_texp = 5*ones(mp.Nsbp,1); % Exposure time for each sub-bandpass (seconds)
    tb.info.sbp_texp_probe = 5*ones(mp.Nsbp,1); % Exposure time for each sub-bandpass when probing (seconds)
end