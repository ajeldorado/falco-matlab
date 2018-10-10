% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute the column of the DM control Jacobian for the
% specified DM and actuator. Uses the compact model.
% 
% Modified on 2018-05-23 by A.J. Riggs to include weights and eliminate the
% zeroing out of the command first (to allow a non-zero starting point).
% Created on 2018-03-28 by A.J. Riggs.

function JacCol = func_validate_Jacobian_with_compact_model_HLC(iact,whichDM,dV,EunpokedVec,mp, DM,  lambda, normFac, Ein, FPM)
          
    act_sens = 1; %--Default unless overwritten
    stepFac = 1; %--Default unless overwritten
    
    flagSkipAct = true;
    
    if(whichDM==1)
        weight = mp.dm_weights(1);
        mp.dm1.V(iact) = mp.dm1.V(iact) + dV; %--Poke one actuator a tiny amount to stay linear
    
        if(any(any(mp.dm1.compact.inf_datacube(:,:,iact))))
            flagSkipAct = false;
        end
        
    elseif(whichDM==2)
        weight = mp.dm_weights(2);
        mp.dm2.V(iact) = mp.dm2.V(iact) + dV; %--Poke one actuator a tiny amount to stay linear
    
        if(any(any(mp.dm2.compact.inf_datacube(:,:,iact))))
            flagSkipAct = false;
        end
        
    elseif(whichDM==8)
        weight = mp.dm_weights(8);
        mp.dm8.V(iact) = mp.dm8.V(iact) + dV; %--Poke one actuator a tiny amount to stay linear
        if(isfield(mp.dm8,'act_sens'))
            act_sens = mp.dm8.act_sens;
        end
        
        if(any(any(mp.dm8.compact.inf_datacube(:,:,iact))))
            flagSkipAct = false;
        end
        
    elseif(whichDM==9)
        
        if(isfield(mp.dm9,'stepFac')==false)
            stepFac = 20;%10; %--Adjust the step size in the Jacobian, then divide back out. Used for helping counteract effect of discretization.
        else
            stepFac = mp.dm9.stepFac;
        end
        
        weight = mp.dm_weights(9);
        mp.dm9.V(iact) = mp.dm9.V(iact) + stepFac*dV; %--Poke one actuator a tiny amount to stay linear
        if(isfield(mp.dm9,'act_sens'))
            act_sens = mp.dm9.act_sens;
        end
        
        if(any(any(mp.dm9.compact.inf_datacube(:,:,iact))))
            flagSkipAct = false;
        end
    end

    
    if(flagSkipAct)
        JacCol = zeros(size(mp.F4.corr.inds)); 
    else
        Epoked = model_compact_HLC(mp, DM,  lambda, normFac, Ein, FPM);
        JacCol = act_sens*weight/stepFac*(Epoked(mp.F4.corr.inds)-EunpokedVec)/dV; %--column of the Jacobian for actuator # iact
    end
end %--END OF FUNCTION
