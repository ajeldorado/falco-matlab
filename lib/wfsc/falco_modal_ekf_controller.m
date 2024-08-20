% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Output control commands from modal controller.
%
% INPUTS
% ------
% mp : structure of model parameters
% cvar : structure of controller variables
%
% OUTPUTS
% -------
% dDM : structure of the delta control commands separated by DM number.
%       Also contains the updated array of tied actuator pairs
% cvar : structure of controller variables

function [dDM, cvar] = falco_modal_ekf_controller(mp, cvar)
    
    cvar = falco_ctrl_setup(mp,cvar);
    duVec = -cvar.u_hat;
    [~,dDM] = falco_ctrl_wrapup(mp,cvar,duVec);
       
end %--END OF FUNCTION
