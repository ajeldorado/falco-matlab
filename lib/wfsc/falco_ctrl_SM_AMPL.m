% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function for regularized linear least-squares control with boundary
% constraints imposed by AMPL and Gurobi.
% -This function performs an empirical grid search over these parameters:
%  a) a scalar coefficient for the regularization matrix
%  b) a scalar gain for the final DM command.
%
% -This code is 
%
%
% REVISION HISTORY:
% -Created on 2018-10-15 by A.J. Riggs.

function [dDM,cvar] = falco_ctrl_SM_AMPL(mp,cvar)



    %% Make vector of variable combinations to try
    vals_list = mp.ctrl.log10regVec; %--dimensions: [1 x length(mp.ctrl.muVec) ]
    %         vals_list = allcomb(mp.ctrl.muVec,mp.ctrl.dmfacVec).'; %--dimensions: [2 x length(mp.ctrl.muVec)*length(mp.ctrl.dmfacVec) ]

    Nvals = max(size(vals_list,2));
    
    %% Temporarily store computed DM commands so that the best one does not have to be re-computed
    
%     if(any(mp.dm_ind==1)); dDM1V_store = zeros(mp.dm1.Nact,mp.dm1.Nact,Nvals); end
%     if(any(mp.dm_ind==2)); dDM2V_store = zeros(mp.dm2.Nact,mp.dm2.Nact,Nvals); end
%     if(any(mp.dm_ind==8)); dDM8V_store = zeros(mp.dm8.NactTotal,Nvals); end
%     if(any(mp.dm_ind==9)); dDM9V_store = zeros(mp.dm9.NactTotal,Nvals); end
    
    %--Storage array for average normalized intensity
    Inorm_list = zeros(Nvals,1);    
    
    %% Make the vector of input, total control commands
    if(any(mp.dm_ind==1));  u1 = mp.dm1.V(mp.dm1.act_ele);  else;  u1 = [];  end
    if(any(mp.dm_ind==2));  u2 = mp.dm2.V(mp.dm2.act_ele);  else;  u2 = [];  end
    if(any(mp.dm_ind==8));  u8 = mp.dm8.V(mp.dm8.act_ele);  else;  u8 = [];  end
    if(any(mp.dm_ind==9));  u9 = mp.dm9.V(mp.dm9.act_ele);  else;  u9 = [];  end
    u = [u1; u2; u8; u9]; %--column vector
    Nele = length(u);
    
    %--Get the indices of each DM's command within the full command
    if(any(mp.dm_ind==1));  u1dummy = 1*ones(mp.dm1.Nele,1);  else;  u1dummy = [];  end
    if(any(mp.dm_ind==2));  u2dummy = 2*ones(mp.dm2.Nele,1);  else;  u2dummy = [];  end
    if(any(mp.dm_ind==8));  u8dummy = 8*ones(mp.dm8.Nele,1);  else;  u8dummy = [];  end
    if(any(mp.dm_ind==9));  u9dummy = 9*ones(mp.dm9.Nele,1);  else;  u9dummy = [];  end
    u_dm_guide = [u1dummy; u2dummy; u8dummy; u9dummy];

    %%
    

    %% Constraints on Actuation
    %--Put constraints in same format (du isolated on one side)
    % du <= du_max
    % du >= -du_max
    % du <= u_ub - u0
    % du >= u_lb - u0
    % % u0+du <= u_ub
    % % u0+du >= u_lb

    %--Constraints: Lower bounds on absolute control commands
    if(any(mp.dm_ind==1));  u1_LB = -1*mp.dm1.maxAbsV*ones(1,mp.dm1.Nele);  else;  u1_LB = [];  end
    if(any(mp.dm_ind==2));  u2_LB = -1*mp.dm2.maxAbsV*ones(1,mp.dm2.Nele);  else;  u2_LB = [];  end
    if(any(mp.dm_ind==8));  u8_LB = mp.dm8.Vmin*ones(1,mp.dm8.Nele);  else;  u8_LB = [];  end
    if(any(mp.dm_ind==9));  u9_LB = mp.dm9.Vmin*ones(1,mp.dm9.Nele);  else;  u9_LB = [];  end
    u_LB = [u1_LB, u2_LB, u8_LB, u9_LB].'; %--column vector
    %--Constraints: Upper bounds on absolute control commands
    if(any(mp.dm_ind==1));  u1_UB = mp.dm1.maxAbsV*ones(1,mp.dm1.Nele);  else;  u1_UB = [];  end
    if(any(mp.dm_ind==2));  u2_UB = mp.dm2.maxAbsV*ones(1,mp.dm2.Nele);  else;  u2_UB = [];  end
    if(any(mp.dm_ind==8));  u8_UB = mp.dm8.Vmax*ones(1,mp.dm8.Nele);  else;  u8_UB = [];  end
    if(any(mp.dm_ind==9));  u9_UB = mp.dm9.Vmax*ones(1,mp.dm9.Nele);  else;  u9_UB = [];  end
    u_UB = [u1_UB, u2_UB, u8_UB, u9_UB].'; %--column vector

    %--Constraints: Lower and upper bounds on delta control commands
    if(any(mp.dm_ind==1));  du1_UB = mp.dm1.maxAbsdV*ones(1,mp.dm1.Nele);  else;  du1_UB = [];  end
    if(any(mp.dm_ind==2));  du2_UB = mp.dm2.maxAbsdV*ones(1,mp.dm2.Nele);  else;  du2_UB = [];  end
    if(any(mp.dm_ind==8));  du8_UB = mp.dm8.maxAbsdV*ones(1,mp.dm8.Nele);  else;  du8_UB = [];  end
    if(any(mp.dm_ind==9));  du9_UB = mp.dm9.maxAbsdV*ones(1,mp.dm9.Nele);  else;  du9_UB = [];  end
    du_UB = [du1_UB, du2_UB, du8_UB, du9_UB].'; %--column vector
    du_LB = -du_UB;

    %--Constraints: Lower and uper bounds on new control commands (including starting point)
    du_LB_total = u_LB - u;
    du_UB_total = u_UB - u;

    %--Combined delta constraints (this is what AMPL sees)
    du_LB_comb = max(du_LB_total,du_LB);
    du_UB_comb = min(du_UB_total,du_UB);
    
    
    %%
    fn_GstarG = ['/home/ajriggs/Repos/falco-matlab/data/Jacobians/' sprintf('GstarG_s%d_t%d.DAT',mp.SeriesNum,mp.TrialNum)];
    %    fn_GstarG = [mp.path.falco 'data/Jacobians/' sprintf('GstarG_s%d_t%d.DAT',mp.SeriesNum,mp.TrialNum)];
% %    % dlmwrite(fn_GstarG,cvar.GstarG_wsum,'delimiter', ' ', 'precision', '%.16g');% ,'precision','%.16g'); %--Write out matrix from Matlab
   
   %--Writing with fprintf is Wayyyyy faster than with dlmwrite.
   fprintf('Saving Jacobian*Jacobian to a file for AMPL to read in...')
   fid = fopen(fn_GstarG,'w');
   fprintf(fid, '%.16g ',cvar.GstarG_wsum);
   fclose(fid);
   fprintf('done. Time = %.2fs\n',toc);
   
   %%
   %keyboard
   
    %% Loop over AMPL calls for different regularizations

%         for ni = 1:Nvals
%             [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_AMPL_func(mp,cvar,vals_list,ni,Nele,du_LB_comb,du_UB_comb,u_dm_guide,fn_GstarG);
%         end %--End of loop over regularizations
    
    if(mp.flagParfor)
        parfor ni = 1:Nvals
            [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_AMPL_func(mp,cvar,vals_list,ni,Nele,du_LB_comb,du_UB_comb,u_dm_guide,fn_GstarG);
        end %--End of loop over regularizations
    else
        for ni = 1:Nvals
            [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_AMPL_func(mp,cvar,vals_list,ni,Nele,du_LB_comb,du_UB_comb,u_dm_guide,fn_GstarG);
        end %--End of loop over regularizations
    end


% % %--Doesn't work because "Using license file "/home/ajriggs/bin/amplapi/matlab/../lib/../../ampl.lic". is not serializable "
%     ampl_shared = falco_ctrl_SM_AMPL_func_a(cvar,Nele,du_LB_comb,du_UB_comb,fn_GstarG);   
%     for ni = 1:Nvals
%        ampl{ni} = ampl_shared; 
%     end
%     
%     if(mp.flagParfor)
%         parfor ni = 1:Nvals
%             [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_AMPL_func_b(mp,vals_list,ni,u_dm_guide,ampl)
%         end %--End of loop over regularizations
%     else
%         for ni = 1:Nvals
%             %[Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_AMPL_func(mp,cvar,vals_list,ni,Nele,du_LB_comb,du_UB_comb,u_dm_guide,fn_GstarG);
%         end %--End of loop over regularizations
%     end
      
    
    

    %% Print out results to the command line if more than one used
    if(numel(mp.ctrl.log10regVec)>1)

        fprintf('\nlog10reg:\t');
        for ni=1:Nvals
            fprintf('%.2e\t',vals_list(1,ni))
        end

        fprintf('\nInorm:  \t')
        for ni=1:Nvals
            fprintf('%.2e\t',Inorm_list(ni))
        end
        fprintf('\n')

    end
        
    

    %% Find the best regularization value based on the best normalized intensity.
    [cvar.cMin,indBest] = min(Inorm_list(:));
    dDM = dDM_cells{indBest};
    
%     if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest); end
%     if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest); end
%     if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V_store(:,indBest); end
%     if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V_store(:,indBest); end
  
    cvar.log10regUsed = vals_list(1,indBest);
    fprintf('Empirically chosen log10reg = %.2f\t   gives %4.2e normalized intensity.\n',cvar.log10regUsed,cvar.cMin)    

    %cp.muBest(Itr) = vals_list(1,indBest);
    %fprintf('Empirically chosen mu = %.2e\t   gives %4.2e contrast.\n',cp.muBest(Itr),cvar.cMin)    


    
    %%
   


end %--END OF FUNCTION
