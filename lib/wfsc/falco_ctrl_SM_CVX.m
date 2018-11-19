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

function [dDM,cvar] = falco_ctrl_SM_CVX(mp,cvar)



    %% Make vector of variable combinations to try
    vals_list = mp.ctrl.log10regVec; %--dimensions: [1 x length(mp.ctrl.muVec) ]
    %         vals_list = allcomb(mp.ctrl.muVec,mp.ctrl.dmfacVec).'; %--dimensions: [2 x length(mp.ctrl.muVec)*length(mp.ctrl.dmfacVec) ]

    Nvals = max(size(vals_list,2));
    
    %--Storage array for average normalized intensity
    Inorm_list = zeros(Nvals,1);   
        
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
    u = [u1; u2; u3; u4; u5; u6; u7; u8; u9]; %--column vector
    NeleAll = length(u);

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
    u_guide = [u1dummy; u2dummy; u3dummy;  u4dummy;  u5dummy; u6dummy;  u7dummy;  u8dummy; u9dummy];

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
    if(any(mp.dm_ind==5));  u5_LB = mp.dm5.Vmin*ones(1,mp.dm5.Nele);  else;  u5_LB = [];  end
    if(any(mp.dm_ind==8));  u8_LB = mp.dm8.Vmin*ones(1,mp.dm8.Nele);  else;  u8_LB = [];  end
    if(any(mp.dm_ind==9));  u9_LB = mp.dm9.Vmin*ones(1,mp.dm9.Nele);  else;  u9_LB = [];  end
    u_LB = [u1_LB, u2_LB, u5_LB, u8_LB, u9_LB].'; %--column vector
    %--Constraints: Upper bounds on absolute control commands
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
    du_LB_total = u_LB - u;
    du_UB_total = u_UB - u;

    %--Combined delta constraints (this is what AMPL sees)
    du_LB_comb = max(du_LB_total,du_LB);
    du_UB_comb = min(du_UB_total,du_UB);
    
    %% Diagonal of the regularization matrix
    
    
    %--The entire regularization matrix will be normalized based on the max
    %response of DMs 1 and 2
    temp = diag(cvar.GstarG_wsum);
    diagNormVal = max(temp(u_guide==1 | u_guide==2));
    
    %--Initialize the regularization matrix
    RegMatDiag = diagNormVal*ones(NeleAll,1);
    
    %--Change regularization values for non-standard DMs
    if(any(mp.dm_ind==8));  RegMatDiag(u_guide==8) = mp.ctrl.relReg8 ;  end
    if(any(mp.dm_ind==9));  RegMatDiag(u_guide==9) = mp.ctrl.relReg9 ;  end
    
    
    %% Loop over CVX calls for different regularizations
    for ni = 1:Nvals
        [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_CVX_func(mp,cvar,vals_list,ni,du_LB_comb,du_UB_comb,u,u_guide,RegMatDiag);
    end %--End of loop over regularizations
    
    

    
    
    

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
  
    cvar.log10regUsed = vals_list(1,indBest);
    fprintf('Empirically chosen log10reg = %.2f\t   gives %4.2e normalized intensity.\n',cvar.log10regUsed,cvar.cMin)    


end %--END OF FUNCTION





%     %% Code for AMPL version
%     fn_GstarG = ['/home/ajriggs/Repos/falco-matlab/data/Jacobians/' sprintf('GstarG_s%d_t%d.DAT',mp.SeriesNum,mp.TrialNum)];
%     %    fn_GstarG = [mp.path.falco 'data/Jacobians/' sprintf('GstarG_s%d_t%d.DAT',mp.SeriesNum,mp.TrialNum)];
% % %    % dlmwrite(fn_GstarG,cvar.GstarG_wsum,'delimiter', ' ', 'precision', '%.16g');% ,'precision','%.16g'); %--Write out matrix from Matlab
%    
%    %--Writing with fprintf is Wayyyyy faster than with dlmwrite.
%    fprintf('Saving Jacobian*Jacobian to a file for AMPL to read in...')
%    fid = fopen(fn_GstarG,'w');
%    fprintf(fid, '%.16g ',cvar.GstarG_wsum);
%    fclose(fid);
%    fprintf('done. Time = %.2fs\n',toc);
%  
%     %% Loop over AMPL calls for different regularizations
% 
% %         for ni = 1:Nvals
% %             [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_AMPL_func(mp,cvar,vals_list,ni,Nele,du_LB_comb,du_UB_comb,u_guide,fn_GstarG);
% %         end %--End of loop over regularizations
%     
%     if(mp.flagParfor)
%         parfor ni = 1:Nvals
%             [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_AMPL_func(mp,cvar,vals_list,ni,NeleAll,du_LB_comb,du_UB_comb,u_guide,fn_GstarG);
%         end %--End of loop over regularizations
%     else
%         for ni = 1:Nvals
%             [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_AMPL_func(mp,cvar,vals_list,ni,NeleAll,du_LB_comb,du_UB_comb,u_guide,fn_GstarG);
%         end %--End of loop over regularizations
%     end
