function ev = falco_est_ekf_maintenance(mp, ev, varargin

%% This stuff has been copy-pasted

Itr = ev.Itr;
whichDM = mp.est.probe.whichDM;

if ~isa(mp.est.probe, 'Probe')
    error('mp.est.probe must be an instance of class Probe')
end


% Augment which DMs are used if the probing DM isn't used for control.
if whichDM == 1 && ~any(mp.dm_ind == 1)
    mp.dm_ind = [mp.dm_ind(:); 1];
elseif whichDM == 2 && ~any(mp.dm_ind == 2)
    mp.dm_ind = [mp.dm_ind(:); 2];
end

%--Select number of actuators across based on chosen DM for the probing
if whichDM == 1
    Nact = mp.dm1.Nact;
elseif whichDM == 2
    Nact = mp.dm2.Nact;
end

%--Store the initial DM commands
% sfr: what is this doing
if any(mp.dm_ind == 1);  DM1Vnom = mp.dm1.V;  else; DM1Vnom = zeros(size(mp.dm1.V)); end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2);  DM2Vnom = mp.dm2.V;  else; DM2Vnom = zeros(size(mp.dm2.V)); end % The 'else' block would mean we're only using DM1

% Initialize output arrays
ev.imageArray = zeros(mp.FendNeta, mp.Fend.Nxi, 1, mp.Nsbp);
ev.Eest = zeros(mp.Fend.corr.Npix, mp.Nsbp*mp.compact.star.count);
ev.IincoEst = zeros(mp.Fend.corr.Npix, mp.Nsbp*mp.compact.star.count);
ev.IprobedMean = 0;
ev.Im = zeros(mp.Fend.Neta, mp.Fend.Nxi);
if whichDM == 1;  ev.dm1.Vall = zeros(mp.dm1.Nact, mp.dm1.Nact, 1+2*Npairs, mp.Nsbp);  end
if whichDM == 2;  ev.dm2.Vall = zeros(mp.dm2.Nact, mp.dm2.Nact, 1+2*Npairs, mp.Nsbp);  end

%% Initialize items
if mp.Itr == 1
    mp = initialize_ekf_maintenance(mp);
end

%% Get Drift Command
mp = falco_drift_injection(mp);

%% Get dither command
% Is this if else necessary?

if any(mp.dm_ind == 1);  DM1Vdither = normrnd(0,ev.dither,[mp.dm1.Nact mp.dm1.Nact]);  else; DM1Vdither = zeros(size(mp.dm1.V)); end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2);  DM2Vdither = normrnd(0,ev.dither,[mp.dm1.Nact mp.dm1.Nact]);  else; DM2Vdither = zeros(size(mp.dm2.V)); end % The 'else' block would mean we're only using DM1


%% Set total command for estimator image
% TODO: need to save these commands for each iteration separately
% TODO: what is controller command?

mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + DM1Vdither + DM1Vdrift + DM1V_controller_result);
mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + DM2Vdither + DM1Vdrift + DM2V_controller_result);

%% Get images


for iSubband = 1:mp.Nsbp
    I0 = falco_get_sbp_image(mp, iSubband) * tb.info.sbp_texp(si) * tb.info.PSFpeaks(si);


end

%% Perform the estimation

for iSubband = 1:mp.Nsbp


end

end