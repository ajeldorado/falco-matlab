% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to handle the output command vectors from the controller.
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - cvar = structure of controller variables
% - duVec = the vector of delta commands computed by the controller
%
% OUTPUTS
% - mp = structure of model parameters
% - dDM = structure of the delta control commands separated by DM number.
%         Also contains the updated array of tied actuator pairs
%
% REVISION HISTORY
% - Created on 2019-02-13 by A.J. Riggs.
% - Modified on 2019-02-25 by A.J. Riggs to save the delta steps.
% - Modified on 2019-03-26 by A.J. Riggs to include tied actuators.
% - Modified on 2019-06-25 by A.J. Riggs to tie actuators that violate the neighbor rule.
% - Modified on 2019-09-26 by A.J. Riggs to handle DM1 and DM2 actuator 
% constraints outside this function in a more user-robust way.

function [mp,dDM] = falco_ctrl_wrapup(mp,cvar,duVec)

%--Initialize delta DM commands
if(any(mp.dm_ind==1)); dDM.dDM1V = zeros(mp.dm1.Nact,mp.dm1.Nact); end
if(any(mp.dm_ind==2)); dDM.dDM2V = zeros(mp.dm2.Nact,mp.dm2.Nact); end
if(any(mp.dm_ind==3)); dDM.dDM3V = zeros(mp.dm3.NactTotal,1); end
if(any(mp.dm_ind==4)); dDM.dDM4V = zeros(mp.dm4.NactTotal,1); end
if(any(mp.dm_ind==5)); dDM.dDM5V = zeros(mp.dm5.NactTotal,1); end
if(any(mp.dm_ind==6)); dDM.dDM6V = zeros(mp.dm6.NactTotal,1); end
if(any(mp.dm_ind==7)); dDM.dDM7V = zeros(mp.dm7.NactTotal,1); end
if(any(mp.dm_ind==8)); dDM.dDM8V = zeros(mp.dm8.NactTotal,1); end
if(any(mp.dm_ind==9)); dDM.dDM9V = zeros(mp.dm9.NactTotal,1); end

%--Parse the delta command vector by DM
if(any(mp.dm_ind==1));  dDM.dDM1V(mp.dm1.act_ele) = mp.dm1.weight*duVec(cvar.uLegend==1);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==2));  dDM.dDM2V(mp.dm2.act_ele) = mp.dm2.weight*duVec(cvar.uLegend==2);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==3));  dDM.dDM3V(mp.dm3.act_ele) = mp.dm3.weight*duVec(cvar.uLegend==3);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==4));  dDM.dDM4V(mp.dm4.act_ele) = mp.dm4.weight*duVec(cvar.uLegend==4);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==5));  dDM.dDM5V(mp.dm5.act_ele) = mp.dm5.weight*duVec(cvar.uLegend==5);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==6));  dDM.dDM6V(mp.dm6.act_ele) = mp.dm6.weight*duVec(cvar.uLegend==6);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==7));  dDM.dDM7V(mp.dm7.act_ele) = mp.dm7.weight*duVec(cvar.uLegend==7);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==8));  dDM.dDM8V(mp.dm8.act_ele) = mp.dm8.weight*duVec(cvar.uLegend==8);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==9));  dDM.dDM9V(mp.dm9.act_ele) = mp.dm9.weight*duVec(cvar.uLegend==9);  end % Parse the command vector to get component for DM and apply the DM's weight

%--Enforce tied actuator pair commands. 
% Assign command of first actuator to the second as well.
if(any(mp.dm_ind==8))
    for ti=1:size(mp.dm8.tied,1)
        dDM.dDM8V(mp.dm8.tied(ti,2)) = dDM.dDM8V(mp.dm8.tied(ti,1));
    end
end
if(any(mp.dm_ind==9))
    for ti=1:size(mp.dm9.tied,1)
        dDM.dDM9V(mp.dm9.tied(ti,2)) = dDM.dDM9V(mp.dm9.tied(ti,1));
    end
end

if(any(mp.dm_ind==1));  mp.dm1.dV = dDM.dDM1V;  end % Store the delta DM command
if(any(mp.dm_ind==2));  mp.dm2.dV = dDM.dDM2V;  end % Store the delta DM command
if(any(mp.dm_ind==3));  mp.dm3.dV = dDM.dDM3V;  end % Store the delta DM command
if(any(mp.dm_ind==4));  mp.dm4.dV = dDM.dDM4V;  end % Store the delta DM command
if(any(mp.dm_ind==5));  mp.dm5.dV = dDM.dDM5V;  end % Store the delta DM command
if(any(mp.dm_ind==6));  mp.dm6.dV = dDM.dDM6V;  end % Store the delta DM command
if(any(mp.dm_ind==7));  mp.dm7.dV = dDM.dDM7V;  end % Store the delta DM command
if(any(mp.dm_ind==8));  mp.dm8.dV = dDM.dDM8V;  end % Store the delta DM command
if(any(mp.dm_ind==9));  mp.dm9.dV = dDM.dDM9V;  end % Store the delta DM command

%--Combine the delta command with the previous command
if(any(mp.dm_ind==1));  mp.dm1.V = cvar.DM1Vnom + dDM.dDM1V;  end
if(any(mp.dm_ind==2));  mp.dm2.V = cvar.DM2Vnom + dDM.dDM2V;  end
if(any(mp.dm_ind==3));  mp.dm3.V = cvar.DM3Vnom + dDM.dDM3V;  end
if(any(mp.dm_ind==4));  mp.dm4.V = cvar.DM4Vnom + dDM.dDM4V;  end
if(any(mp.dm_ind==5));  mp.dm5.V = cvar.DM5Vnom + dDM.dDM5V;  end
if(any(mp.dm_ind==6));  mp.dm6.V = cvar.DM6Vnom + dDM.dDM6V;  end
if(any(mp.dm_ind==7));  mp.dm7.V = cvar.DM7Vnom + dDM.dDM7V;  end
if(any(mp.dm_ind==8));  mp.dm8.V = cvar.DM8Vnom + dDM.dDM8V;  end
if(any(mp.dm_ind==9));  mp.dm9.V = cvar.DM9Vnom + dDM.dDM9V;  end

end %--END OF FUNCTION