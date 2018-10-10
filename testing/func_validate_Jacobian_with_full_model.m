% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute the column of the DM control Jacobian for the
% specified DM and actuator. Uses the full model.
% 
% Modified on 2018-05-23 by A.J. Riggs to include weights and eliminate the
% zeroing out of the command first (to allow a non-zero starting point).
% Created on 2018-03-28 by A.J. Riggs.


function JacCol = func_validate_Jacobian_with_full_model(iact,whichDM,dV,EunpokedVec,mp, DM, modvar)
         
    if(whichDM==1)
        weight = mp.dm_weights(1);
        mp.dm1.V(iact) = mp.dm1.V(iact) + dV; %--Poke one actuator a tiny amount to stay linear
    elseif(whichDM==2)
        weight = mp.dm_weights(2);
        mp.dm2.V(iact) = mp.dm2.V(iact) + dV; %--Poke one actuator a tiny amount to stay linear
    elseif(whichDM==9)
        weight = mp.dm_weights(9);
        mp.dm9.V(iact) = mp.dm9.V(iact) + dV; %--Poke one actuator a tiny amount to stay linear
    end
    Epoked = model_full(mp, DM, modvar);

    JacCol = weight*(Epoked(mp.F4.corr.inds)-EunpokedVec)/dV; %--column of the Jacobian for actuator # iact

end %--END OF FUNCTION
