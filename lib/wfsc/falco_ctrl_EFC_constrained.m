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

function dDM = falco_ctrl_EFC_constrained(DM, cvar)

%% Optimization using CVX

if( any(DM.dm_ind==1) && any(DM.dm_ind==2) && length(DM.dm_ind)==2 ) %--DMs 1 and 2 only
    DM1VnomVec = DM.dm1.V(DM.dm1.act_ele);
    DM2VnomVec = DM.dm2.V(DM.dm2.act_ele);
    
    cvx_begin quiet
        % cvx_precision low
        cvx_solver Mosek
        variables maxContrast u1(DM.dm1.Nele) u2(DM.dm2.Nele)
        minimize (maxContrast)
        subject to
            ([u1; u2].' * cvar.GstarG_wsum + cvar.RealGstarEab_wsum.') * [u1; u2] <= maxContrast

            u1 <= DM.dm1.dVpvMax/2
            u1 >= -DM.dm1.dVpvMax/2
            u2 <= DM.dm2.dVpvMax/2
            u2 >= -DM.dm2.dVpvMax/2
            
            u1 + DM1VnomVec <= DM.dm1.maxAbsV
            u1 + DM1VnomVec >= -DM.dm1.maxAbsV
            u2 + DM2VnomVec <= DM.dm2.maxAbsV
            u2 + DM2VnomVec >= -DM.dm2.maxAbsV
            
    cvx_end
end


%% Reformat the solution for the output

if(any(DM.dm_ind==1))
    dDM1V = zeros(DM.dm1.Nact);
    dDM1V(DM.dm1.act_ele) = u1;
    dDM.dDM1V = dDM1V;
end

if(any(DM.dm_ind==2))
    dDM2V = zeros(DM.dm2.Nact, DM.dm2.Nact);
    dDM2V(DM.dm2.act_ele) = u2;
    dDM.dDM2V = dDM2V;
end

if(any(DM.dm_ind==3))
    dDM3V = zeros(DM.dm3.Nact);
    dDM3V(DM.dm3.act_ele) = u3;
    dDM.dDM3V = dDM3V;
end




end