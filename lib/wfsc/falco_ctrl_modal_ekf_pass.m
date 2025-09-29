% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Peform EFC according to a pre-defined (planned) scheme.
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

function [dDM, cvar] = falco_ctrl_modal_ekf_pass(mp, cvar)
    % Set gain to 0 for first N iterations
    if cvar.Itr<mp.ctrl.start_iteration; gain = 0; else gain = 1; end
    cvar.du_hat = gain*cvar.du_hat;
    cvar = falco_ctrl_setup(mp, cvar);

    duVec = -mean(cvar.du_hat,2);
    
    [mp,dDM] = falco_ctrl_wrapup(mp,cvar,duVec);

    Itotal = falco_get_expected_summed_image(mp,cvar);
    InormMean = mean(Itotal(mp.Fend.corr.maskBool));
    dDM.Itotal = Itotal;

end %--END OF FUNCTION