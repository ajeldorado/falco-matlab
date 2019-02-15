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
        
    %%
    %--Initialize delta DM commands
    %--Parse the command vector by DM
    %--Combine the delta command with the previous command
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
    
    
    %% Loop over CVX calls for different regularizations
    for ni = 1:Nvals
        [Inorm_list(ni),dDM_cells{ni}] = falco_ctrl_SM_CVX_func(mp,cvar,vals_list,ni);
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
