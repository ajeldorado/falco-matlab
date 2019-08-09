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
% - Modified on 2019-04-23 by A.J. Riggs to have an option for a
%   model-based grid search.
% - Modified on 2019-02-13 by A.J. Riggs to use falco_ctrl_setup.m and
%   falco_ctrl_wrapup.m to save a bunch of space.
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

function [InormAvg,thput,dDM] = falco_ctrl_EFC_base(ni,vals_list,nj,valsOmega_list,nk,vals_list_dm9,mp,cvar)


%% Initializations
% Itr = cvar.Itr ;
log10reg = vals_list(1,ni); %--Lagrange multiplier
dmfac = vals_list(2,ni); %--Scaling factor for entire DM command
log10regdm9 = vals_list_dm9(nk);
mp.aux.omega = valsOmega_list(nj);
%--Save starting point for each delta command to be added to.
%--Get the indices of each DM's command vector within the single concatenated command vector
cvar = falco_ctrl_setup(mp,cvar);

%% Define the diagonal of the regularization matrix differently
%--Diagonal of the Weighted Regularization Matrix
% EyeGstarGdiag = [];
% maxDiagGstarG = max(diag(cvar.GstarG_wsum));
% for idm=1:numel(mp.dm_ind)
%     dm_index = mp.dm_ind(idm);
%     dm_weight = 1; %mp.dm_weights(dm_index);
%     if(any(mp.dm_ind==9)); dm_weight = dm9regfac*dm_weight; end
%     EyeGstarGdiag = [EyeGstarGdiag; maxDiagGstarG*dm_weight*ones(cvar.NeleVec(idm),1)];
% end
%% Define the regularization for DM9 differently
reg_diag = cvar.EyeGstarGdiag;
maxDiagGstarG = max(diag(cvar.GstarG_wsum));
reg_diag = reg_diag*10^(log10reg);
ind_dm9 = find(cvar.uLegend==9);
reg_diag(ind_dm9) = 10^(log10regdm9); 

%% Least-squares solution:
% duVec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;
% duVec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\(cvar.RealGstarEab_wsum+0.0*10^(log10reg)*cvar.EyeGstarGdiag.*...
%     [mp.dm1.V(mp.dm1.act_ele);mp.dm2.V(mp.dm2.act_ele)]);
if(any(mp.dm_ind==9))
    
    if(any(mp.aux.dm9OnlyItr_arr==cvar.Itr))
        aux_ind_mat = zeros(size(cvar.GstarG_wsum));
        sz_dm9 = numel(mp.dm9.act_ele);
        cvar.EyeGstarGdiag = cvar.EyeGstarGdiag(end-sz_dm9+1:end);
        cvar.GstarG_wsum = cvar.GstarG_wsum(end-sz_dm9+1:end,end-sz_dm9+1:end);
        cvar.RealGstarEab_wsum = cvar.RealGstarEab_wsum(end-sz_dm9+1:end);
        reg_diag = reg_diag(end-sz_dm9+1:end);
        vec_dm_ele = [mp.dm9.V(mp.dm9.act_ele)];
    else
        vec_dm_ele = [mp.dm1.V(mp.dm1.act_ele);mp.dm2.V(mp.dm2.act_ele);mp.dm9.V(mp.dm9.act_ele)];
    end
else
    vec_dm_ele = [mp.dm1.V(mp.dm1.act_ele);mp.dm2.V(mp.dm2.act_ele)];
end
if mp.aux.flagOmega==1 && cvar.Itr>=mp.aux.firstOmegaItr
    duVec = -dmfac*(diag(reg_diag) + cvar.GstarG_wsum - 10^mp.aux.omega * cvar.GcptransGcp_wsum)...
        \(cvar.RealGstarEab_wsum...
        +mp.aux.gamma*10^(log10reg)*cvar.EyeGstarGdiag.*...
        vec_dm_ele);
else
    duVec = -dmfac*(diag(reg_diag) + mp.aux.gamma*cvar.EyeGstarGdiag + cvar.GstarG_wsum)...
        \(cvar.RealGstarEab_wsum+mp.aux.gamma*cvar.EyeGstarGdiag.*...
        vec_dm_ele);
end

if(any(mp.aux.dm9OnlyItr_arr==cvar.Itr) && any(mp.dm_ind==9)) 
    duVec = [zeros(numel(mp.dm1.act_ele)+numel(mp.dm2.act_ele),1);duVec];
end

% dDMvec = -dmfac*(diag(EyeGstarGdiag)/mu + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;
% % dDMvec = -dmfac*(diag(cvar.EyeGstarGdiag)/mu + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;


%% Parse the command vector by DM and assign the output commands
[mp,dDM] = falco_ctrl_wrapup(mp,cvar,duVec);

%% Take images and compute average intensity in dark hole

if(mp.ctrl.flagUseModel) %--Perform a model-based grid search using the compact model
    if(mp.flagFiber)
        %--Not available yet
    else
        Itotal = falco_get_expected_summed_image(mp,cvar);
        InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
    end
else %--Perform an empirical grid search with actual images from the testbed or full model
    if(mp.flagFiber)
        IfiberTotal = falco_get_summed_image_fiber(mp);
        InormAvg = mean(max(max(IfiberTotal)));
    else
        Itotal = falco_get_summed_image(mp);
        InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
        [mp,thput] = falco_compute_thput(mp);
    end
end
        

end %--END OF FUNCTION





