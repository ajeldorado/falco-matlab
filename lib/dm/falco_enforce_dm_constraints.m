% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Apply various constraints on the actuators of the deformable 
% mirrors. This is for conventional DMs, DMs 1 and 2.
%
% Constraints:
% 1) Min and max bounds
% 2) Pinned/railed actuators
% 3) Neighbor rule and tied actuators
%
% INPUTS
% ------
% dm : structure of DM parameters
%
% OUTPUTS
% -------
% dm : structure of DM parameters

function dm = falco_enforce_dm_constraints(dm)

% 1) Find actuators that below min value or above max value.
% Any actuators reaching those limits are added to the pinned actuator list.

% Min voltage limit
Vtotal = dm.V+dm.biasMap;
indVoltageTooLow = find(Vtotal < dm.Vmin);
dm.pinned = [dm.pinned; indVoltageTooLow]; % augment the column vector of pinned actuators' linear indices
dm.Vpinned = [dm.Vpinned; (dm.Vmin*ones(size(indVoltageTooLow))-dm.biasMap(indVoltageTooLow))];

% Max voltage limit
indVoltageTooHigh = find(Vtotal > dm.Vmax);
dm.pinned = [dm.pinned; indVoltageTooHigh]; % augment the column vector of pinned actuators' linear indices
dm.Vpinned = [dm.Vpinned; (dm.Vmax*ones(size(indVoltageTooHigh))-dm.biasMap(indVoltageTooHigh))];

% 2) Enforce bounds at pinned actuators
dm.V(dm.pinned) = dm.Vpinned; 

% 3) Enforce neighbor rule and tied actuators at same time (actually
% iterated between the two).
tieMat = falco_convert_dm_tie_pairs_into_matrix(dm.tied, dm.Nact);
vlat = dm.dVnbr;
vdiag = dm.dVnbr;
maxiter = 1000;
vquant = 0; % LSB in volts
Vtotal = dm.V + dm.biasMap;
Vtotal = ConstrainDM.constrain_dm(Vtotal, dm.biasMap, tieMat, dm.Vmax, vlat, vdiag, vquant, maxiter);
dm.V = Vtotal - dm.biasMap;

% % 3) Find which actuators violate the DM neighbor rule. (This restricts 
% % the maximum voltage between an actuator and each of its 8 neighbors.) 
% % Add those actuator pairs to the list of tied actuators.
% if dm.flagNbrRule
%     [dm.V, indPair1] = falco_dm_neighbor_rule(dm.V, dm.dVnbr, dm.Nact);
%     dm.tied = [dm.tied; indPair1]; %--Tie together actuators violating the neighbor rule
% end
%     
% % 4) Enforce tied actuator pairs
% % In each pair of tied actuators, assign the command for the 1st actuator to that of the 2nd actuator
% if ~isempty(dm.tied)
%     dm.V(dm.tied(:, 2)) = dm.V(dm.tied(:, 1));% + dm.biasMap(dm.tied(:, 1)) - dm.biasMap(dm.tied(:, 2));
% end

end
