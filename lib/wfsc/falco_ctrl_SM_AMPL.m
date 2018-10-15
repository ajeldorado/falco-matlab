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
    vals_list = mp.ctrl.muVec; %--dimensions: [1 x length(mp.ctrl.muVec) ]
    %         vals_list = allcomb(mp.ctrl.muVec,mp.ctrl.dmfacVec).'; %--dimensions: [2 x length(mp.ctrl.muVec)*length(mp.ctrl.dmfacVec) ]

    Nvals = max(size(vals_list,2));
    
    %% Temporarily store computed DM commands so that the best one does not have to be re-computed
    if(any(mp.dm_ind==1)); dDM1V_store = zeros(mp.dm1.Nact,mp.dm1.Nact,Nvals); end
    if(any(mp.dm_ind==2)); dDM2V_store = zeros(mp.dm2.Nact,mp.dm2.Nact,Nvals); end
    if(any(mp.dm_ind==8)); dDM8V_store = zeros(mp.dm8.NactTotal,Nvals); end
    if(any(mp.dm_ind==9)); dDM9V_store = zeros(mp.dm9.NactTotal,Nvals); end
    
    %--Storage array for average normalized intensity
    Inorm_list = zeros(Nvals,1);    
    
    %% Make the vector of input, total control commands
    if(any(mp.dm_ind==1));  u1 = mp.dm1.V(mp.dm1.act_ele);  else;  u1 = [];  end
    if(any(mp.dm_ind==2));  u2 = mp.dm2.V(mp.dm2.act_ele);  else;  u2 = [];  end
    if(any(mp.dm_ind==8));  u8 = mp.dm8.V(mp.dm8.act_ele);  else;  u8 = [];  end
    if(any(mp.dm_ind==9));  u9 = mp.dm9.V(mp.dm9.act_ele);  else;  u9 = [];  end
    u = [u1, u2, u8, u9].'; %--column vector
    Nele = length(u);
    
    %--Get the indices of each DM's command within the full command
    if(any(mp.dm_ind==1));  u1dummy = 1*ones(mp.dm1.Nele,1);  else;  u1dummy = [];  end
    if(any(mp.dm_ind==2));  u2dummy = 2*ones(mp.dm2.Nele,1);  else;  u2dummy = [];  end
    if(any(mp.dm_ind==8));  u8dummy = 8*ones(mp.dm8.Nele,1);  else;  u8dummy = [];  end
    if(any(mp.dm_ind==9));  u9dummy = 9*ones(mp.dm9.Nele,1);  else;  u9dummy = [];  end
    u_dm_guide = [u1dummy, u2dummy, u8dummy, u9dummy];

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
    if(any(mp.dm_ind==1));  u1_LB = -1*mp.dm1.maxAbsV*ones(1,mp.dm1.Nele);  else;  du1_LB = [];  end
    if(any(mp.dm_ind==2));  u2_LB = -1*mp.dm2.maxAbsV*ones(1,mp.dm2.Nele);  else;  du2_LB = [];  end
    if(any(mp.dm_ind==8));  u8_LB = mp.dm8.Vmin*ones(1,mp.dm8.Nele);  else;  du8_LB = [];  end
    if(any(mp.dm_ind==9));  u9_LB = mp.dm9.Vmin*ones(1,mp.dm9.Nele);  else;  du9_LB = [];  end
    u_LB = [u1_LB, u2_LB, u8_LB, u9_LB].'; %--column vector
    %--Constraints: Upper bounds on absolute control commands
    if(any(mp.dm_ind==1));  u1_UB = mp.dm1.maxAbsV*ones(1,mp.dm1.Nele);  else;  du1_UB = [];  end
    if(any(mp.dm_ind==2));  u2_UB = mp.dm2.maxAbsV*ones(1,mp.dm2.Nele);  else;  du2_UB = [];  end
    if(any(mp.dm_ind==8));  u8_UB = mp.dm8.Vmax*ones(1,mp.dm8.Nele);  else;  du8_UB = [];  end
    if(any(mp.dm_ind==9));  u9_UB = mp.dm9.Vmax*ones(1,mp.dm9.Nele);  else;  du9_UB = [];  end
    u_UB = [u1_UB, u2_UB, u8_UB, u9_UB].'; %--column vector

    %--Constraints: Lower and upper bounds on delta control commands
    if(any(mp.dm_ind==1));  du1_UB = mp.dm1.maxAbsdV*ones(1,mp.dm1.Nele);  else;  du1_LB = [];  end
    if(any(mp.dm_ind==2));  du2_UB = mp.dm2.maxAbsdV*ones(1,mp.dm2.Nele);  else;  du2_LB = [];  end
    if(any(mp.dm_ind==8));  du8_UB = mp.dm8.maxAbsdV*ones(1,mp.dm8.Nele);  else;  du8_LB = [];  end
    if(any(mp.dm_ind==9));  du9_UB = mp.dm9.maxAbsdV*ones(1,mp.dm9.Nele);  else;  du9_LB = [];  end
    du_UB = [du1_UB, du2_UB, du8_UB, du9_UB].'; %--column vector
    du_LB = -du_UB;

    %--Constraints: Lower and uper bounds on new control commands (including starting point)
    du_LB_total = u_LB - u;
    du_UB_total = u_UB - u;

    %--Combined delta constraints (this is what AMPL sees)
    du_LB_comb = max(du_LB_total,du_LB);
    du_UB_comb = min(du_UB_total,du_UB);

    
    
    %% Loop over AMPL calls
    
    
    for ni = 1:Nvals
        
         %% Interface with AMPL

        % Create an AMPL instance
        ampl = AMPL;
        % Display version
        ampl.eval('option version;')

        ampl.eval( sprintf('param Nele := %d;',Nele) );
        ampl.eval('set Acts := setof {q in 1..Nele by 1} q;');                                                        
        ampl.eval(sprintf('param log10reg := %.3f;',log10reg));

        %--Transfer lower bound on delta command
        ampl.eval('param ampl_du_LB {q in Acts};');
        ampl_du_LB = ampl.getParameter('ampl_du_LB');
        ampl_du_LB.setValues(du_LB_comb); %--Transfer values from Matlab to AMPL

        %--Transfer upper bound on delta command
        ampl.eval('param ampl_du_UB {q in Acts};');
        ampl_du_UB = ampl.getParameter('ampl_du_UB');
        ampl_du_UB.setValues(du_UB_comb); %--Transfer values from Matlab to AMPL

        %--Define the delta control vector with bounds
        ampl.eval('var du {q in Acts};'); 
        ampl.eval('var du {q in Acts} >=ampl_du_LB , <=ampl_du_UB, := 0;');


        %--Control Matrices
        ampl.eval('param GstarG {x in Acts,y in Acts};');
        GstarG = ampl.getParameter('GstarG');
        GstarG.setValues(cvar.GstarG_wsum); %--Transfer values from Matlab to AMPL

        ampl.eval('param RealGstarEab {q in Acts};');
        RealGstarEab = ampl.getParameter('RealGstarEab');
        RealGstarEab.setValues(cvar.RealGstarEab_wsum); %--Transfer values from Matlab to AMPL

        %--Regularization Matrix: Read in the diagonal, make the matrix of zeros, and fill in the diagonal.
        ampl.eval('param EyeGstarGdiag {q in Acts};');
        EyeGstarGdiag = ampl.getParameter('EyeGstarGdiag');
        EyeGstarGdiag.setValues(cvar.EyeGstarGdiag); %--Transfer values from Matlab to AMPL
        ampl.eval('param RegMat {x in Acts, y in Acts};');
        ampl.eval('let {x in Acts, y in Acts} RegMat[x,y] := 0;');
        ampl.eval('let {x in Acts, y in Acts: x==y} RegMat[x,y] := EyeGstarGdiag[x];');

        ampl.eval('param QuadMat {x in Acts, y in Acts};');
        ampl.eval('let {x in Acts, y in Acts} QuadMat[x,y] := RegMat[x,y] + GstarG[x,y];');


        ampl.eval('var costStep1 {q in Acts};'); 
        ampl.eval('var J;');  %--The cost function
        ampl.eval('minimize total_cost: J; ')     

        ampl.eval('subject to costStep1_def {x in Acts, y in Acts}:  costStep1[y] = sum {x in Acts} du[x]*(GstarG[x,y] + 10^(log10reg)*RegMat[x,y]);');
        ampl.eval('subject to costStep2_def {q in Acts}:  J = sum {q in Acts}  (costStep1[q] + RealGstarEab[q])*du[q];')

        ampl.eval('option solver gurobi;');
        ampl.eval('option gurobi_options "outlev=1 presolve=0 lpmethod=2 barhomogeneous=1 crossover=0 threads=28";');

        % ([u1; u2].' * (cvar.GstarG_wsum + 10^log10reg*maxDiagGstarG*eye(size(cvar.GstarG_wsum))) + cvar.RealGstarEab_wsum.') * [u1; u2] <= maxContrast
        % dDMvec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;

        %--Retrieve results
        totalcost = ampl.getObjective('total_cost');
        fprintf('Objective is: %f\n', totalcost.value);

        % Get the values of the variable Buy in a dataframe object
        du_temp = ampl.getVariable('du');
        du_out = du_temp.getValues;
        % Print them
        % du_out

        %--Close AMPL
        ampl.reset();
        
        
        %--Assign the output commands
        if(any(mp.dm_ind==1)); dDM1V = zeros(mp.dm1.Nact,mp.dm1.Nact);  dDM1V(mp.dm1.act_ele) = du_out(u_dm_guide==1);  dDM1V_store(:,:,ni) = dDM1V; end
        if(any(mp.dm_ind==2)); dDM2V = zeros(mp.dm2.Nact,mp.dm2.Nact);  dDM2V(mp.dm2.act_ele) = du_out(u_dm_guide==2);  dDM2V_store(:,:,ni) = dDM2V; end
        if(any(mp.dm_ind==8)); dDM8V_store(u_dm_guide==8,ni) = dDM8V; end
        if(any(mp.dm_ind==9)); dDM9V_store(u_dm_guide==9,ni) = dDM9V; end
        
        %% Take images to empirically check contrast at that setting
        Itotal = falco_get_summed_image(mp);
        Inorm_list(ni) = mean(Itotal(mp.F4.corr.maskBool));
     

    
    end
    

    %--Print out results to the command line if more than one used
    if(numel(mp.ctrl.log10regVec)>1)

        fprintf('\nlog10reg:\t\t');
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

    if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest); end
    if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest); end
    if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V_store(:,indBest); end
    if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V_store(:,indBest); end
  
    log10regBest = vals_list(1,indBest);
    fprintf('Empirically chosen log10reg = %.2e\t   gives %4.2e normalized intensity.\n',log10regBest,cvar.cMin)    

    %cp.muBest(Itr) = vals_list(1,indBest);
    %fprintf('Empirically chosen mu = %.2e\t   gives %4.2e contrast.\n',cp.muBest(Itr),cvar.cMin)    

    
    
    
    %%
    %%
    %%
    
%%
    
    
%     % Temporarily store computed DM commands so that the best one does not have to be re-computed
%     if(any(mp.dm_ind==1)); dDM1V_store = zeros(mp.dm1.Nact,mp.dm1.Nact,Nvals); end
%     if(any(mp.dm_ind==2)); dDM2V_store = zeros(mp.dm2.Nact,mp.dm2.Nact,Nvals); end
%     if(any(mp.dm_ind==8)); dDM8V_store = zeros(mp.dm8.NactTotal,Nvals); end
%     if(any(mp.dm_ind==9)); dDM9V_store = zeros(mp.dm9.NactTotal,Nvals); end
% 
% 
% %     if(mp.flagParfor) %--Parallelized %--DOES NOT WORK WITH CVX. THE CVX PROGRAMS OVERWRITE EACH OTHER.
% %         parfor ni = 1:Nvals
% %             [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_constrained_func(ni,vals_list,DM, cvar);%(ni,vals_list,DM, cp, cvar, mp);
% %             if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
% %             if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
% %             if(any(mp.dm_ind==9)); dDM9V_store(:,ni) = dDM_temp.dDM9V; end
% %         end
% %     else %--Not Parallelized
%         for ni = 1:Nvals
%             [Inorm_list(ni),dDM_temp] = falco_ctrl_SM_AMPL_func(ni,vals_list, cvar, mp);
%             if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
%             if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
%             if(any(mp.dm_ind==8)); dDM8V_store(:,ni) = dDM_temp.dDM8V; end
%             if(any(mp.dm_ind==9)); dDM9V_store(:,ni) = dDM_temp.dDM9V; end
%         end
% %     end
% 
% 
%     %--Print out results to the command line if more than one used
% if(numel(mp.ctrl.muVec)>1)
% 
%     fprintf('\nMu:\t\t');
%     for ni=1:Nvals
%         fprintf('%.2e\t',vals_list(1,ni))
%     end
% 
%     fprintf('\nInorm:  \t')
%     for ni=1:Nvals
%         fprintf('%.2e\t',Inorm_list(ni))
%     end
%     fprintf('\n')
% 
% end



%     %--Find the best Lagrange multiplier value based on the best contrast.
%     [cvar.cMin,indBest] = min(Inorm_list(:));
% 
%     if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest); end
%     if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest); end
%     if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V_store(:,indBest); end
%     if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V_store(:,indBest); end
%   
%     muBest = vals_list(1,indBest);
%     fprintf('Empirically chosen mu = %.2e\t   gives %4.2e contrast.\n',muBest,cvar.cMin)    
% 
%     %cp.muBest(Itr) = vals_list(1,indBest);
%     %fprintf('Empirically chosen mu = %.2e\t   gives %4.2e contrast.\n',cp.muBest(Itr),cvar.cMin)    


end %--END OF FUNCTION
