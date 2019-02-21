% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function for regularized linear least-squares control (EFC).
% -This function computes the DM commands and new Inorm for one set of these parameters:
%  a) a scalar coefficient for the regularization matrix
%  b) a scalar gain for the final DM command.
%
% -This code is based on electric field conjugation (EFC) as described 
% by Give'on et al. SPIE 2011.
%
%
% REVISION HISTORY: 
% -Modified on 2018-11-12 by A.J. Riggs to clean up code, especially to 
%   remove the large commented-out blocks.
% -Modified on 2018-07-24 by A.J. Riggs to switch from the Lagrange multiplier to the Tikhonov regularization.
% -Modified on 2018-02-06 by A.J. Riggs to be parallelized with parfor.
%   Called by a higher function. 
% -Modified by A.J. Riggs on October 11, 2017 to allow easier mixing of
%   which DMs are used and to also do a grid search over the gain of the 
%   overall DM command. 
% -Modified from hcil_ctrl_checkMuEmp.m by A.J. Riggs on August 31, 2016
% -Created at Princeton on 19 Feb 2015 by A.J. Riggs

%--Return values:
%  Measured average normalized intensity
%  DM commands
function [InormAvg,dDM] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar)

%% Initializations
log10reg = vals_list(1,ni); %--Lagrange multiplier
dmfac = vals_list(2,ni); %--Scaling factor for entire DM command

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
cvar.u_guide = [u1dummy; u2dummy; u3dummy;  u4dummy;  u5dummy; u6dummy;  u7dummy;  u8dummy; u9dummy];

%% Least-squares solution:
duVec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;

%% Parse the command vector by DM

%--Initialize delta DM commands
if(any(mp.dm_ind==1)); dDM.dDM1V = zeros(mp.dm1.Nact,mp.dm1.Nact); end
if(any(mp.dm_ind==2)); dDM.dDM2V = zeros(mp.dm2.Nact,mp.dm2.Nact); end
% if(any(mp.dm_ind==1)); dDM.dDM1V = zeros(mp.dm1.NactTotal,1); end
% if(any(mp.dm_ind==2)); dDM.dDM2V = zeros(mp.dm2.NactTotal,1); end
if(any(mp.dm_ind==3)); dDM.dDM3V = zeros(mp.dm3.NactTotal,1); end
if(any(mp.dm_ind==4)); dDM.dDM4V = zeros(mp.dm4.NactTotal,1); end
if(any(mp.dm_ind==5)); dDM.dDM5V = zeros(mp.dm5.NactTotal,1); end
if(any(mp.dm_ind==6)); dDM.dDM6V = zeros(mp.dm6.NactTotal,1); end
if(any(mp.dm_ind==7)); dDM.dDM7V = zeros(mp.dm7.NactTotal,1); end
if(any(mp.dm_ind==8)); dDM.dDM8V = zeros(mp.dm8.NactTotal,1); end
if(any(mp.dm_ind==9)); dDM.dDM9V = zeros(mp.dm9.NactTotal,1); end

%--Parse the command vector by DM
if(any(mp.dm_ind==1));  dDM.dDM1V(mp.dm1.act_ele) = mp.dm_weights(1)*duVec(cvar.u_guide==1);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==2));  dDM.dDM2V(mp.dm2.act_ele) = mp.dm_weights(2)*duVec(cvar.u_guide==2);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==3));  dDM.dDM3V(mp.dm3.act_ele) = mp.dm_weights(3)*duVec(cvar.u_guide==3);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==4));  dDM.dDM4V(mp.dm4.act_ele) = mp.dm_weights(4)*duVec(cvar.u_guide==4);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==5));  dDM.dDM5V(mp.dm5.act_ele) = mp.dm_weights(5)*duVec(cvar.u_guide==5);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==6));  dDM.dDM6V(mp.dm6.act_ele) = mp.dm_weights(6)*duVec(cvar.u_guide==6);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==7));  dDM.dDM7V(mp.dm7.act_ele) = mp.dm_weights(7)*duVec(cvar.u_guide==7);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==8));  dDM.dDM8V(mp.dm8.act_ele) = mp.dm_weights(8)*duVec(cvar.u_guide==8);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==9));  dDM.dDM9V(mp.dm9.act_ele) = mp.dm_weights(9)*duVec(cvar.u_guide==9);  end % Parse the command vector to get component for DM and apply the DM's weight

%% Combine the delta command with the previous command

if(any(mp.dm_ind==1));  mp.dm1.V = cvar.DM1Vnom + dDM.dDM1V;  end
if(any(mp.dm_ind==2));  mp.dm2.V = cvar.DM2Vnom + dDM.dDM2V;  end
if(any(mp.dm_ind==3));  mp.dm3.V = cvar.DM3Vnom + dDM.dDM3V;  end
if(any(mp.dm_ind==4));  mp.dm4.V = cvar.DM4Vnom + dDM.dDM4V;  end
if(any(mp.dm_ind==5));  mp.dm5.V = cvar.DM5Vnom + dDM.dDM5V;  end
if(any(mp.dm_ind==6));  mp.dm6.V = cvar.DM6Vnom + dDM.dDM6V;  end
if(any(mp.dm_ind==7));  mp.dm7.V = cvar.DM7Vnom + dDM.dDM7V;  end
if(any(mp.dm_ind==8));  mp.dm8.V = cvar.DM8Vnom + dDM.dDM8V;  end
if(any(mp.dm_ind==9));  mp.dm9.V = cvar.DM9Vnom + dDM.dDM9V;  end

%% Take images and compute average intensity in dark hole

[Itotal, IfiberTotal] = falco_get_summed_image(mp);
if(mp.flagFiber)
    InormAvg = mean(IfiberTotal);
else
    InormAvg = mean(Itotal(mp.F4.corr.maskBool));
end
    
end %--END OF FUNCTION