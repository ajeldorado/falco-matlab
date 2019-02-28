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
% - Modified on 2019-02-13 by A.J. Riggs to use falco_ctrl_setup.m and
% falco_ctrl_wrapup.m to save a bunch of space.
% - Modified on 2018-11-12 by A.J. Riggs to clean up code, especially to 
%   remove the large commented-out blocks.
% - Modified on 2018-07-24 by A.J. Riggs to switch from the Lagrange multiplier to the Tikhonov regularization.
% - Modified on 2018-02-06 by A.J. Riggs to be parallelized with parfor.
%   Called by a higher function. 
% - Modified by A.J. Riggs on October 11, 2017 to allow easier mixing of
%   which DMs are used and to also do a grid search over the gain of the 
%   overall DM command. 
% - Modified from hcil_ctrl_checkMuEmp.m by A.J. Riggs on August 31, 2016
% - Created at Princeton on 19 Feb 2015 by A.J. Riggs


%--Return values:
%  Measured average normalized intensity
%  DM commands

function [InormAvg,dDM] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar)


%% Initializations
% Itr = cvar.Itr ;
log10reg = vals_list(1,ni); %--Lagrange multiplier
dmfac = vals_list(2,ni); %--Scaling factor for entire DM command

%--Save starting point for each delta command to be added to.
%--Get the indices of each DM's command vector within the single concatenated command vector
cvar = falco_ctrl_setup(mp,cvar);

%% Define the diagonal of the regularization matrix differently
% % %--Diagonal of the Weighted Regularization Matrix
% % EyeGstarGdiag = [];
% % maxDiagGstarG = max(diag(cvar.GstarG_wsum));
% % for idm=1:numel(mp.dm_ind)
% %     dm_index = mp.dm_ind(idm);
% %     dm_weight = 1; %mp.dm_weights(dm_index);
% %     if(any(mp.dm_ind==9)); dm_weight = dm9regfac*dm_weight; end
% %     EyeGstarGdiag = [EyeGstarGdiag; maxDiagGstarG*dm_weight*ones(cvar.NeleVec(idm),1)];
% % end

%% Least-squares solution:
% duVec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;
duVec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\(cvar.RealGstarEab_wsum+0.0*10^(log10reg)*cvar.EyeGstarGdiag.*...
    [mp.dm1.V(mp.dm1.act_ele);mp.dm2.V(mp.dm2.act_ele)]);
% dDMvec = -dmfac*(diag(EyeGstarGdiag)/mu + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;
% % dDMvec = -dmfac*(diag(cvar.EyeGstarGdiag)/mu + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;


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

%% Parse the command vector by DM and assign the output commands
[mp,dDM] = falco_ctrl_wrapup(mp,cvar,duVec);

%% Take images and compute average intensity in dark hole
Itotal = falco_get_summed_image(mp);
InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
        

end %--END OF FUNCTION





