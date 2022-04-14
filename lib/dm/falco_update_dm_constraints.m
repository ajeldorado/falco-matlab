% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Update the lists of which DM actuators are pinned or tied.
%
% dm.tied stays the same because those are electrically tied actuators.
% dm.comovingGroups is used instead for new tied actuators from the
% neighbor rule because those behave differently (they have a constant,
% nonzero voltage offset) and should not be used as inputs to
% ConstrainDM.constrain_dm().

function dm = falco_update_dm_constraints(dm)

    % Get 2-D map of electrically tied actuators, as specified by dm.tied.
    % (Electrically tied here means that they always share the same voltage.)
    constraintMap = falco_convert_dm_tie_pairs_into_matrix(dm.tied, dm.Nact);
    
    % Add in known dead actuators.
    % Railed ones will be handled by ActLimits.maplimits()
    constraintMap(dm.dead) = -1;

    % Update constraint map based on current voltages
    % Actuators become frozen (pinned) if they are at the upper or lower
    % bounds. Actuators become tied (at a specific offset) if they hit the
    % neighbor rule limit.
    Vtotal = dm.V + dm.biasMap;
    [dm.pinned, dm.comovingGroups] = ActLimits.maplimits(Vtotal, dm, constraintMap);
    dm.Vpinned = dm.V(dm.pinned);

end