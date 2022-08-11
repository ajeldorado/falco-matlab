
function [mp,ev] = falco_drift_injection(mp,ev)
    % TODO: eventually move drift to different function and put in main loop
    % before estimator
switch lower(mp.drift.type)
    case{'rand_walk'}
        if any(mp.dm_drift_ind == 1)
            mp.dm1.V_drift = mp.dm1.V_drift + normrnd(0,mp.drift.magnitude,[mp.dm1.Nact mp.dm1.Nact]);  
%             mp.dm1.V = mp.dm1.V + mp.dm1.V_drift;
        else
            mp.dm1.V_drift = zeros(size(mp.dm1.V));
        end % The 'else' block would mean we're only using DM2
        
        if any(mp.dm_drift_ind == 2)  
            mp.dm2.V_drift = mp.dm2.V_drift + normrnd(0,mp.drift.magnitude,[mp.dm1.Nact mp.dm1.Nact]);  
%             mp.dm2.V = mp.dm2.V + mp.dm2.V_drift;
        else 
            mp.dm2.V_drift = zeros(size(mp.dm2.V)); 
        end % The 'else' block would mean we're only using DM1

end
if any(mp.est.itr_reset==ev.Itr) == true
    if ~isfield(mp.dm1,'dV'); mp.dm1.dV = zeros(mp.dm1.Nact);end
    if ~isfield(mp.dm2,'dV'); mp.dm2.dV = zeros(mp.dm2.Nact);end
    efc_command = get_dm_command_vector(mp,mp.dm1.dV, mp.dm2.dV);
    
    for iSubband = 1:1:mp.Nsbp
        ev.x_hat(:,iSubband) = ev.x_hat(:,iSubband) + (ev.G_tot(:,:,iSubband)*ev.e_scaling(iSubband))*sqrt(mp.tb.info.sbp_texp(iSubband))*efc_command;
    end
    mp.dm1.V_shift = mp.dm1.dV;
    mp.dm2.V_shift = mp.dm2.dV;

%     mp.dm1.V_dz = mp.dm1.V_dz + mp.dm1.dV;
%     mp.dm2.V_dz = mp.dm2.V_dz + mp.dm2.dV;

    mp.dm1.dV = zeros(size(mp.dm1.V_dz));
    mp.dm2.dV = zeros(size(mp.dm2.V_dz));
end
    % TODO: save each drift command
end



% TODO: this should not be duplicated

function comm_vector = get_dm_command_vector(mp,command1, command2)

if any(mp.dm_ind == 1); comm1 = command1(mp.dm1.act_ele);  else; comm1 = []; end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2); comm2 = command2(mp.dm2.act_ele);  else; comm2 = []; end
comm_vector = [comm1;comm2];

end


