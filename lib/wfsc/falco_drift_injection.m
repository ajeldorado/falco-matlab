function mp = falco_drift_injection(mp)
    % TODO: eventually move drift to different function and put in main loop
    % before estimator
    if any(mp.dm_drift_ind == 1);  mp.dm1.V_drift = mp.dm1.V_drift + normrnd(0,mp.drift.magnitude,[mp.dm1.Nact mp.dm1.Nact]);  else; mp.dm1.V_drift = zeros(size(mp.dm1.V)); end % The 'else' block would mean we're only using DM2
    if any(mp.dm_drift_ind == 2);  mp.dm2.V_drift = mp.dm2.V_drift + normrnd(0,mp.drift.magnitude,[mp.dm1.Nact mp.dm1.Nact]);  else; mp.dm2.V_drift = zeros(size(mp.dm2.V)); end % The 'else' block would mean we're only using DM1
    
    % TODO: save each drift command
end



