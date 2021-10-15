% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% 
% Set the DM voltage and then apply constraints.
%
% INPUTS
% ------
% dm : structure of DM parameters
% Vnew : new voltage command to assign to dm.V
%
% OUTPUTS
% -------
% dm : structure of DM parameters

function dm = falco_set_constrained_voltage(dm, Vnew)

    dm.V = Vnew;
    dm = falco_enforce_dm_constraints(dm);

end