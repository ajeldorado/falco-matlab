 function [InormAvg,dDM] = falco_ctrl_SM_CVX_func(mp,cvar,vals_list,ni)
%  function [InormAvg,dDM] = falco_ctrl_SM_CVX_func(mp,cvar,vals_list,ni,du_LB_comb,du_UB_comb,RegMatDiag)
    
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
        
%         if(any(mp.dm_ind==1)); dDM.dDM1V = zeros(mp.dm1.Nact,mp.dm1.Nact);  dDM.dDM1V(mp.dm1.act_ele) = mp.dm1.weight*duVec(uLegend==1);   end
%         if(any(mp.dm_ind==2)); dDM.dDM2V = zeros(mp.dm2.Nact,mp.dm2.Nact);  dDM.dDM2V(mp.dm2.act_ele) = mp.dm2.weight*duVec(uLegend==2);   end
%         if(any(mp.dm_ind==8)); dDM.dDM8V = zeros(mp.dm8.NactTotal,1);  dDM.dDM8V(mp.dm8.act_ele) = mp.dm8.weight*duVec(uLegend==8);  end
%         if(any(mp.dm_ind==9)); dDM.dDM9V = zeros(mp.dm9.NactTotal,1);  dDM.dDM9V(mp.dm9.act_ele) = mp.dm9.weight*duVec(uLegend==9);  end        
%         
%         %--Update the DM commands by adding the new control signal
%         if(any(mp.dm_ind==1))
%             mp.dm1.dV = dDM.dDM1V;
%             mp.dm1.V = mp.dm1.V + mp.dm1.dV; 
%         end
%         if(any(mp.dm_ind==2))
%             mp.dm2.dV = dDM.dDM2V;
%             mp.dm2.V = mp.dm2.V + mp.dm2.dV; 
%         end
%         if(any(mp.dm_ind==8))
%             mp.dm8.dV = dDM.dDM8V(:);
%             mp.dm8.V = mp.dm8.V + mp.dm8.dV(:);
%         end
%         if(any(mp.dm_ind==9))
%             mp.dm9.dV = dDM.dDM9V;
%             mp.dm9.V(:) = mp.dm9.V(:) + mp.dm9.dV(:);
%         end

        
        %% Take images to empirically check contrast at that setting
        Itotal = falco_get_summed_image(mp);
        InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
        
    end %--END OF NESTED FUNCTION
