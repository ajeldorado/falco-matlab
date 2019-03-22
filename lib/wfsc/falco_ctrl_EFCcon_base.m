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
% -Modified on 2019-03-18 by A.J. Riggs to be for constrained EFC.
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

function [InormAvg,dDM] = falco_ctrl_EFCcon_base(ni,vals_list,mp,cvar)


%% Initializations
% Itr = cvar.Itr ;
log10reg = vals_list(1,ni); %--Lagrange multiplier
dmfac = vals_list(2,ni); %--Scaling factor for entire DM command

%--Save starting point for each delta command to be added to.
%--Get the indices of each DM's command vector within the single concatenated command vector
cvar = falco_ctrl_setup(mp,cvar);

%% Constraints on Actuation
%--Put constraints in same format (du isolated on one side)
% du <= du_max
% du >= -du_max
% du <= u_ub - u0
% du >= u_lb - u0
% % u0+du <= u_ub
% % u0+du >= u_lb

%--Constraints: Lower bounds on total control commands
if(any(mp.dm_ind==1));  u1_LB = -1*mp.dm1.maxAbsV*ones(1,mp.dm1.Nele);  else;  u1_LB = [];  end
if(any(mp.dm_ind==2));  u2_LB = -1*mp.dm2.maxAbsV*ones(1,mp.dm2.Nele);  else;  u2_LB = [];  end
if(any(mp.dm_ind==5));  u5_LB = mp.dm5.Vmin*ones(1,mp.dm5.Nele);  else;  u5_LB = [];  end
if(any(mp.dm_ind==8));  u8_LB = mp.dm8.Vmin*ones(1,mp.dm8.Nele);  else;  u8_LB = [];  end
if(any(mp.dm_ind==9));  u9_LB = mp.dm9.Vmin*ones(1,mp.dm9.Nele);  else;  u9_LB = [];  end
u_LB = [u1_LB, u2_LB, u5_LB, u8_LB, u9_LB].'; %--column vector
%--Constraints: Upper bounds on total control commands
if(any(mp.dm_ind==1));  u1_UB = mp.dm1.maxAbsV*ones(1,mp.dm1.Nele);  else;  u1_UB = [];  end
if(any(mp.dm_ind==2));  u2_UB = mp.dm2.maxAbsV*ones(1,mp.dm2.Nele);  else;  u2_UB = [];  end
if(any(mp.dm_ind==5));  u5_UB = mp.dm5.Vmax*ones(1,mp.dm5.Nele);  else;  u5_UB = [];  end
if(any(mp.dm_ind==8));  u8_UB = mp.dm8.Vmax*ones(1,mp.dm8.Nele);  else;  u8_UB = [];  end
if(any(mp.dm_ind==9));  u9_UB = mp.dm9.Vmax*ones(1,mp.dm9.Nele);  else;  u9_UB = [];  end
u_UB = [u1_UB, u2_UB, u5_UB, u8_UB, u9_UB].'; %--column vector

%--Constraints: Lower and upper bounds on delta control commands
if(any(mp.dm_ind==1));  du1_UB = mp.dm1.maxAbsdV*ones(1,mp.dm1.Nele);  else;  du1_UB = [];  end
if(any(mp.dm_ind==2));  du2_UB = mp.dm2.maxAbsdV*ones(1,mp.dm2.Nele);  else;  du2_UB = [];  end
if(any(mp.dm_ind==5));  du5_UB = mp.dm5.maxAbsdV*ones(1,mp.dm5.Nele);  else;  du5_UB = [];  end
if(any(mp.dm_ind==8));  du8_UB = mp.dm8.maxAbsdV*ones(1,mp.dm8.Nele);  else;  du8_UB = [];  end
if(any(mp.dm_ind==9));  du9_UB = mp.dm9.maxAbsdV*ones(1,mp.dm9.Nele);  else;  du9_UB = [];  end
du_UB = [du1_UB, du2_UB, du5_UB, du8_UB, du9_UB].'; %--column vector
du_LB = -du_UB;

%--Constraints: Lower and uper bounds on new control commands (including starting point)
du_LB_total = u_LB - cvar.uVec;
du_UB_total = u_UB - cvar.uVec;

%--Combined delta constraints (this is what AMPL sees)
cvar.du_LB_comb = max(du_LB_total,du_LB);
cvar.du_UB_comb = min(du_UB_total,du_UB);
    
%% Diagonal of the regularization matrix


%--The entire regularization matrix will be normalized based on the max
%response of DMs 1 and 2
temp = diag(cvar.GstarG_wsum);
if(any(cvar.uLegend==1) || any(cvar.uLegend==2))
    diagNormVal = max(temp(cvar.uLegend==1 | cvar.uLegend==2));
else
    diagNormVal = max(temp);
end
%--Initialize the regularization matrix
cvar.RegMatDiag = diagNormVal*ones(cvar.NeleAll,1);

%--NOT YET IMPLEMENTED TO DO ANYTHING
%--Change regularization values for non-standard DMs
mp.ctrl.relReg8 = 1;
mp.ctrl.relReg9 = 1;
if(any(mp.dm_ind==8));  cvar.RegMatDiag(cvar.uLegend==8) = mp.ctrl.relReg8*cvar.RegMatDiag(cvar.uLegend==8) ;  end
if(any(mp.dm_ind==9));  cvar.RegMatDiag(cvar.uLegend==9) = mp.ctrl.relReg9*cvar.RegMatDiag(cvar.uLegend==9) ;  end
    


%% Least-squares solution with regularization:
% duVec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;

 %% Perform constrained optimization with CVX
log10reg = vals_list(1,ni);
%NeleAll = length(uNom);

cvx_begin quiet
%         cvx_precision default
%         cvx_solver Mosek
    variables maxContrast duVec(cvar.NeleAll,1)
    % variables maxContrast u1(mp.dm1.Nele) u2(mp.dm2.Nele) u8(mp.dm8.Nele) u9(mp.dm9.Nele) 
    minimize (maxContrast)
    subject to
        (duVec.' * (cvar.GstarG_wsum + 10.^(log10reg)*diag(cvar.RegMatDiag))  + 2*cvar.RealGstarEab_wsum.') * duVec <= maxContrast

        duVec <= cvar.du_UB_comb
        duVec >= cvar.du_LB_comb

cvx_end

%% Parse the command vector by DM and assign the output commands
[mp,dDM] = falco_ctrl_wrapup(mp,cvar,duVec);

%% Take images and compute average intensity in dark hole
Itotal = falco_get_summed_image(mp);
InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
        

end %--END OF FUNCTION





