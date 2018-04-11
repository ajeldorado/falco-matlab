% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute the column of the DM control Jacobian for the
%specified DM and actuator. Uses the full model.
% 
% Created by A.J. Riggs on 2018-03-28.


function JacCol = func_validate_Jacobian_with_full_model(iact,whichDM,Vfrac,EunpokedVec,mp, DM, modvar)
         
    if(whichDM==1)
        DM.dm1.V = 0*DM.dm1.V; %--Reset the voltage map to zero
        DM.dm1.V(iact) = Vfrac; %--Poke one actuator a tiny amount to stay linear
    elseif(whichDM==2)
        DM.dm2.V = 0*DM.dm2.V; %--Reset the voltage map to zero
        DM.dm2.V(iact) = Vfrac; %--Poke one actuator a tiny amount to stay linear
    end

    Epoked = model_full(mp, DM, modvar);

    JacCol = (Epoked(mp.F4.compact.corr.inds)-EunpokedVec)/Vfrac; %--column of the Jacobian for actuator # iact

end %--END OF FUNCTION














