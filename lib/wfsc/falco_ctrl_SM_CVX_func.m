 function [InormAvg,dDM] = falco_ctrl_SM_CVX_func(mp,cvar,vals_list,ni)
%  function [InormAvg,dDM] = falco_ctrl_SM_CVX_func(mp,cvar,vals_list,ni,du_LB_comb,du_UB_comb,RegMatDiag)
    
 %% Perform constrained optimization with CVX
    log10reg = vals_list(1,ni);
    
    cvx_begin quiet
        variables maxContrast duVec(cvar.NeleAll,1)
        minimize (maxContrast)
        subject to
            (duVec.' * (cvar.GstarG_wsum + 10.^(log10reg)*diag(cvar.RegMatDiag))  + 2*cvar.RealGstarEab_wsum.') * duVec <= maxContrast
            
            duVec <= cvar.du_UB_comb
            duVec >= cvar.du_LB_comb           
    cvx_end

    %% Parse the command vector by DM and assign the output commands
    [mp,dDM] = falco_ctrl_wrapup(mp,cvar,duVec);

        %% Take images to empirically check contrast at that setting
        Itotal = falco_get_summed_image(mp);
        InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
        
    end %--END OF NESTED FUNCTION