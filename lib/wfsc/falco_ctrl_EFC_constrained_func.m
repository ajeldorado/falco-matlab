% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function for constrained EFC control for multiple DMs in broadband.
% Here, CVX solves the quadratic cost function directly with constraints, whereas
% regular EFC would perform a left pseudo-inverse with regularization.
%
% Modified on 2018-02-26 by A.J. Riggs at JPL to change the variable names
% for FALCO.
% Developed on Feb. 21, 2018 by He Sun at Princeton University.

% function dDM = falco_ctrl_EFC_constrained_func(ni, vals_list, DM, cvar)
function [InormAvg,dDM] = falco_ctrl_EFC_constrained_func(ni, vals_list, cvar, mp)

%% Optimization using CVX

maxDiagGstarG = max(diag(cvar.GstarG_wsum));

mu = vals_list(1,ni); %--Lagrange multiplier
% if(any(mp.dm_ind==9))
%     dm9regfac = vals_list(2,ni); %--Scaling factor for DM9
% else
%     dm9regfac = 1;
% end


if( any(mp.dm_ind==1) && any(mp.dm_ind==2) && any(mp.dm_ind==8) && any(mp.dm_ind==9) && length(mp.dm_ind)==4 ) %--DMs 1, 2, 8, 9
    DM1VnomVec = mp.dm1.V(mp.dm1.act_ele);
    DM2VnomVec = mp.dm2.V(mp.dm2.act_ele);
    DM8VnomVec = mp.dm8.V(mp.dm8.act_ele);
    DM9VnomVec = mp.dm9.V(mp.dm9.act_ele);
    
    cvx_begin quiet
        %cvx_precision low
        cvx_solver Mosek
        variables maxContrast u1(mp.dm1.Nele) u2(mp.dm2.Nele) u8(mp.dm8.Nele) u9(mp.dm9.Nele) 
        minimize (maxContrast)
        subject to
            ([u1; u2; u8; u9].' * (cvar.GstarG_wsum)  + cvar.RealGstarEab_wsum.') * [u1; u2; u8; u9] <= maxContrast

            u1 <= mp.dm1.dVpvMax/2
            u1 >= -mp.dm1.dVpvMax/2
            u2 <= mp.dm2.dVpvMax/2
            u2 >= -mp.dm2.dVpvMax/2
            u8 <= mp.dm8.dVpvMax/2
            u8 >= -mp.dm8.dVpvMax/2
            u9 <= mp.dm9.dVpvMax/2
            u9 >= -mp.dm9.dVpvMax/2
            
            u1 + DM1VnomVec <= mp.dm1.maxAbsV
            u1 + DM1VnomVec >= -mp.dm1.maxAbsV
            u2 + DM2VnomVec <= mp.dm2.maxAbsV
            u2 + DM2VnomVec >= -mp.dm2.maxAbsV
            u8 + DM8VnomVec <= mp.dm8.Vmax
            u8 + DM8VnomVec >= mp.dm8.Vmin
            u9 + DM9VnomVec <= mp.dm9.Vmax
            u9 + DM9VnomVec >= mp.dm9.Vmin
    cvx_end

elseif( any(mp.dm_ind==1) && any(mp.dm_ind==2) && any(mp.dm_ind==8) && length(mp.dm_ind)==3 ) %--DMs 1, 2, 8 only
    DM1VnomVec = mp.dm1.V(mp.dm1.act_ele);
    DM2VnomVec = mp.dm2.V(mp.dm2.act_ele);
    DM8VnomVec = mp.dm8.V(mp.dm8.act_ele);
    
    cvx_begin quiet
        %cvx_precision low
        cvx_solver Mosek
        variables maxContrast u1(mp.dm1.Nele) u2(mp.dm2.Nele) u8(mp.dm8.Nele) 
        minimize (maxContrast)
        subject to
            ([u1; u2; u8].' * (cvar.GstarG_wsum)  + cvar.RealGstarEab_wsum.') * [u1; u2; u8] <= maxContrast

            u1 <= mp.dm1.dVpvMax/2
            u1 >= -mp.dm1.dVpvMax/2
            u2 <= mp.dm2.dVpvMax/2
            u2 >= -mp.dm2.dVpvMax/2
            u8 <= mp.dm8.dVpvMax/2
            u8 >= -mp.dm8.dVpvMax/2
            
            
            u1 + DM1VnomVec <= mp.dm1.maxAbsV
            u1 + DM1VnomVec >= -mp.dm1.maxAbsV
            u2 + DM2VnomVec <= mp.dm2.maxAbsV
            u2 + DM2VnomVec >= -mp.dm2.maxAbsV
            u8 + DM8VnomVec <= mp.dm8.Vmax
            u8 + DM8VnomVec >= mp.dm8.Vmin
    cvx_end

elseif( any(mp.dm_ind==1) && any(mp.dm_ind==2) && any(mp.dm_ind==9) && length(mp.dm_ind)==3 ) %--DMs 1, 2, 9 only
    DM1VnomVec = mp.dm1.V(mp.dm1.act_ele);
    DM2VnomVec = mp.dm2.V(mp.dm2.act_ele);
    DM9VnomVec = mp.dm9.V(mp.dm9.act_ele);
    
    cvx_begin quiet
        %cvx_precision low
        cvx_solver Mosek
        variables maxContrast u1(mp.dm1.Nele) u2(mp.dm2.Nele) u9(mp.dm9.Nele) 
        minimize (maxContrast)
        subject to
            ([u1; u2; u9].' * (cvar.GstarG_wsum)  + cvar.RealGstarEab_wsum.') * [u1; u2; u9] <= maxContrast

            u1 <= mp.dm1.dVpvMax/2
            u1 >= -mp.dm1.dVpvMax/2
            u2 <= mp.dm2.dVpvMax/2
            u2 >= -mp.dm2.dVpvMax/2
            u9 <= mp.dm9.dVpvMax/2
            u9 >= -mp.dm9.dVpvMax/2
            
            u1 + DM1VnomVec <= mp.dm1.maxAbsV
            u1 + DM1VnomVec >= -mp.dm1.maxAbsV
            u2 + DM2VnomVec <= mp.dm2.maxAbsV
            u2 + DM2VnomVec >= -mp.dm2.maxAbsV
            u9 + DM9VnomVec <= mp.dm9.Vmax
            u9 + DM9VnomVec >= mp.dm9.Vmin
    cvx_end

elseif( any(mp.dm_ind==1) && any(mp.dm_ind==2) && length(mp.dm_ind)==2 ) %--DMs 1 and 2 only
    DM1VnomVec = mp.dm1.V(mp.dm1.act_ele);
    DM2VnomVec = mp.dm2.V(mp.dm2.act_ele);
    
    cvx_begin quiet
        %cvx_precision low
        cvx_solver Mosek
        variables maxContrast u1(mp.dm1.Nele) u2(mp.dm2.Nele)
        minimize (maxContrast)
        subject to
            ([u1; u2].' * (cvar.GstarG_wsum + 1/mu*maxDiagGstarG*eye(size(cvar.GstarG_wsum))) + cvar.RealGstarEab_wsum.') * [u1; u2] <= maxContrast

            u1 <= mp.dm1.dVpvMax/2
            u1 >= -mp.dm1.dVpvMax/2
            u2 <= mp.dm2.dVpvMax/2
            u2 >= -mp.dm2.dVpvMax/2
            
            u1 + DM1VnomVec <= mp.dm1.maxAbsV
            u1 + DM1VnomVec >= -mp.dm1.maxAbsV
            u2 + DM2VnomVec <= mp.dm2.maxAbsV
            u2 + DM2VnomVec >= -mp.dm2.maxAbsV
            
    cvx_end
end


%% Reformat the solution for the output

if(any(mp.dm_ind==1))
    dDM1V = zeros(mp.dm1.Nact);
    dDM1V(mp.dm1.act_ele) = u1;
    dDM.dDM1V = dDM1V;
end

if(any(mp.dm_ind==2))
    dDM2V = zeros(mp.dm2.Nact, mp.dm2.Nact);
    dDM2V(mp.dm2.act_ele) = u2;
    dDM.dDM2V = dDM2V;
end

if(any(mp.dm_ind==3)) 
    dDM3V = zeros(mp.dm3.Nact);
    dDM3V(mp.dm3.act_ele) = u3;
    dDM.dDM3V = dDM3V;
end

if(any(mp.dm_ind==8))
    dDM8V = zeros(mp.dm8.NactTotal,1);
    dDM8V(mp.dm8.act_ele) = u8;
    dDM.dDM8V = dDM8V;
end

if(any(mp.dm_ind==9))
    dDM9V = zeros(mp.dm9.NactTotal,1);
    dDM9V(mp.dm9.act_ele) = u9;
    dDM.dDM9V = dDM9V;
end


%% Take images to empirically check contrast at that mu value

if(any(mp.dm_ind==1)); mp.dm1.V = mp.dm1.V + dDM1V; end
if(any(mp.dm_ind==2)); mp.dm2.V = mp.dm2.V + dDM2V; end
if(any(mp.dm_ind==8)); mp.dm8.V = mp.dm8.V + dDM8V; end
if(any(mp.dm_ind==9)); mp.dm9.V = mp.dm9.V + dDM9V; end

if(numel(mp.ctrl.muVec)>1)
    Itotal = falco_get_summed_image(mp);
    InormAvg = mean(Itotal(mp.F4.corr.maskBool));
else %--Don't waste time evaluating if only one value is used
    InormAvg = 1; %--dummy variable for output
end
