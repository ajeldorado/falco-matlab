% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Incorporate delta DM commands into the total DM commands.
%
% INPUTS
% ------
% mp : structure of model parameters
% dDM : structure of delta DM commands
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_ctrl_update_dm_commands(mp, dDM)

    %--Update the DM commands by adding the delta control signal
    if(any(mp.dm_ind==1));  mp.dm1.V = mp.dm1.V + dDM.dDM1V;  end
    if(any(mp.dm_ind==2));  mp.dm2.V = mp.dm2.V + dDM.dDM2V;  end
    if(any(mp.dm_ind==3));  mp.dm3.V = mp.dm3.V + dDM.dDM3V;  end
    if(any(mp.dm_ind==4));  mp.dm4.V = mp.dm4.V + dDM.dDM4V;  end
    if(any(mp.dm_ind==5));  mp.dm5.V = mp.dm5.V + dDM.dDM5V;  end
    if(any(mp.dm_ind==6));  mp.dm6.V = mp.dm6.V + dDM.dDM6V;  end
    if(any(mp.dm_ind==7));  mp.dm7.V = mp.dm7.V + dDM.dDM7V;  end
    if(any(mp.dm_ind==8));  mp.dm8.V = mp.dm8.V + dDM.dDM8V;  end
    if(any(mp.dm_ind==9));  mp.dm9.V = mp.dm9.V + dDM.dDM9V;  end

    %%--Save the delta from the previous command
    if(any(mp.dm_ind==1));  mp.dm1.dV = dDM.dDM1V;  end
    if(any(mp.dm_ind==2));  mp.dm2.dV = dDM.dDM2V;  end
    if(any(mp.dm_ind==3));  mp.dm3.dV = dDM.dDM3V;  end
    if(any(mp.dm_ind==4));  mp.dm4.dV = dDM.dDM4V;  end
    if(any(mp.dm_ind==5));  mp.dm5.dV = dDM.dDM5V;  end
    if(any(mp.dm_ind==6));  mp.dm6.dV = dDM.dDM6V;  end
    if(any(mp.dm_ind==7));  mp.dm7.dV = dDM.dDM7V;  end
    if(any(mp.dm_ind==8));  mp.dm8.dV = dDM.dDM8V;  end
    if(any(mp.dm_ind==9));  mp.dm9.dV = dDM.dDM9V;  end

end
