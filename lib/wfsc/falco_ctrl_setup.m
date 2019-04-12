% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to set up control variables as vectors for all controllers.
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - cvar = structure of controller variables
%
% OUTPUTS
% - cvar = structure of controller variables
%
% REVISION HISTORY
% - Created on 2019-02-13 by A.J. Riggs.

function cvar = falco_ctrl_setup(mp,cvar)

%--Save starting point for each delta command to be added to.
if(any(mp.dm_ind==1)); cvar.DM1Vnom = mp.dm1.V; end
if(any(mp.dm_ind==2)); cvar.DM2Vnom = mp.dm2.V; end
if(any(mp.dm_ind==3)); cvar.DM3Vnom = mp.dm3.V; end
if(any(mp.dm_ind==4)); cvar.DM4Vnom = mp.dm4.V; end
if(any(mp.dm_ind==5)); cvar.DM5Vnom = mp.dm5.V; end
if(any(mp.dm_ind==6)); cvar.DM6Vnom = mp.dm6.V; end
if(any(mp.dm_ind==7)); cvar.DM7Vnom = mp.dm7.V; end
if(any(mp.dm_ind==8)); cvar.DM8Vnom = mp.dm8.V(:); end
if(any(mp.dm_ind==9)); cvar.DM9Vnom = mp.dm9.V; end

%% Make the vector of input, total control commands
if(any(mp.dm_ind==1));  u1 = mp.dm1.V(mp.dm1.act_ele);  else;  u1 = [];  end
if(any(mp.dm_ind==2));  u2 = mp.dm2.V(mp.dm2.act_ele);  else;  u2 = [];  end
if(any(mp.dm_ind==3));  u3 = mp.dm3.V(mp.dm3.act_ele);  else;  u3 = [];  end
if(any(mp.dm_ind==4));  u4 = mp.dm4.V(mp.dm4.act_ele);  else;  u4 = [];  end
if(any(mp.dm_ind==5));  u5 = mp.dm5.V(mp.dm5.act_ele);  else;  u5 = [];  end
if(any(mp.dm_ind==6));  u6 = mp.dm6.V(mp.dm6.act_ele);  else;  u6 = [];  end
if(any(mp.dm_ind==7));  u7 = mp.dm7.V(mp.dm7.act_ele);  else;  u7 = [];  end
if(any(mp.dm_ind==8));  u8 = mp.dm8.V(mp.dm8.act_ele);  else;  u8 = [];  end
if(any(mp.dm_ind==9));  u9 = mp.dm9.V(mp.dm9.act_ele);  else;  u9 = [];  end
cvar.uVec = [u1; u2; u3; u4; u5; u6; u7; u8; u9]; %--column vector
cvar.NeleAll = length(cvar.uVec);

%--Get the indices of each DM's command within the full command
if(any(mp.dm_ind==1));  u1dummy = 1*ones(mp.dm1.Nele,1);  else;  u1dummy = [];  end
if(any(mp.dm_ind==2));  u2dummy = 2*ones(mp.dm2.Nele,1);  else;  u2dummy = [];  end
if(any(mp.dm_ind==3));  u3dummy = 3*ones(mp.dm3.Nele,1);  else;  u3dummy = [];  end
if(any(mp.dm_ind==4));  u4dummy = 4*ones(mp.dm4.Nele,1);  else;  u4dummy = [];  end
if(any(mp.dm_ind==5));  u5dummy = 5*ones(mp.dm5.Nele,1);  else;  u5dummy = [];  end
if(any(mp.dm_ind==6));  u6dummy = 6*ones(mp.dm6.Nele,1);  else;  u6dummy = [];  end
if(any(mp.dm_ind==7));  u7dummy = 7*ones(mp.dm7.Nele,1);  else;  u7dummy = [];  end
if(any(mp.dm_ind==8));  u8dummy = 8*ones(mp.dm8.Nele,1);  else;  u8dummy = [];  end
if(any(mp.dm_ind==9));  u9dummy = 9*ones(mp.dm9.Nele,1);  else;  u9dummy = [];  end
cvar.uLegend = [u1dummy; u2dummy; u3dummy;  u4dummy;  u5dummy; u6dummy;  u7dummy;  u8dummy; u9dummy];

%--Save starting point for each delta command to be added to.
if(any(mp.dm_ind==1)); cvar.DM1Vnom = mp.dm1.V; end
if(any(mp.dm_ind==2)); cvar.DM2Vnom = mp.dm2.V; end
if(any(mp.dm_ind==3)); cvar.DM3Vnom = mp.dm3.V; end
if(any(mp.dm_ind==4)); cvar.DM4Vnom = mp.dm4.V; end
if(any(mp.dm_ind==5)); cvar.DM5Vnom = mp.dm5.V; end
if(any(mp.dm_ind==6)); cvar.DM6Vnom = mp.dm6.V; end
if(any(mp.dm_ind==7)); cvar.DM7Vnom = mp.dm7.V; end
if(any(mp.dm_ind==8)); cvar.DM8Vnom = mp.dm8.V(:); end
if(any(mp.dm_ind==9)); cvar.DM9Vnom = mp.dm9.V; end

end %--END OF FUNCTION