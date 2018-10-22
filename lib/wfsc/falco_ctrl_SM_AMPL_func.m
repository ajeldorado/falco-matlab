 function [InormAvg,dDM] = falco_ctrl_SM_AMPL_func(mp,cvar,vals_list,ni,Nele,du_LB_comb,du_UB_comb,u_dm_guide,fn_GstarG)
        
        log10reg = vals_list(1,ni);

    
         %%--Interface with AMPL

        % Create an AMPL instance
        setUp();
        ampl{ni} = AMPL;
        % Display version
        %ampl{ni}.eval('option version;')

        ampl{ni}.eval( sprintf('param Nele := %d;',Nele) );
        ampl{ni}.eval('set Acts := setof {q in 1..Nele by 1} q;');                                                        
        ampl{ni}.eval(sprintf('param log10reg := %.3f;',log10reg));

        %--Transfer lower bound on delta command
        ampl{ni}.eval('param ampl_du_LB {q in Acts};');
        ampl_du_LB = ampl{ni}.getParameter('ampl_du_LB');
        ampl_du_LB.setValues(du_LB_comb); %--Transfer values from Matlab to AMPL

        %--Transfer upper bound on delta command
        ampl{ni}.eval('param ampl_du_UB {q in Acts};');
        ampl_du_UB = ampl{ni}.getParameter('ampl_du_UB');
        %ampl_du_UB.setValues(du_UB_comb); %--Transfer values from Matlab to AMPL
        ampl_du_UB.setValues((1:Nele).',du_UB_comb); %--Transfer values from Matlab to AMPL

%          % Get the values of the variable Buy in a dataframe object
%         out_thing = ampl{ni}.getParameter('ampl_du_UB');
%         out_print = out_thing.getValues();
%         % Print them
%         a = out_print.val; whos a
%         b = cell2mat(a); whos b
        
        %--Define the delta control vector with bounds
        ampl{ni}.eval('var du {q in Acts} >=ampl_du_LB[q] , <=ampl_du_UB[q], := 1;');
        %ampl{ni}.eval('var du {q in Acts} >=ampl_du_LB[q] , <=ampl_du_UB[q], := 0;');

        %--Control Matrices
        ampl{ni}.eval('param GstarG {x in Acts,y in Acts};');
        %GstarG = ampl{ni}.getParameter('GstarG');
        %GstarG.setValues(cvar.GstarG_wsum); %--Transfer values from Matlab to AMPL
        ampl{ni}.eval(['read {y in Acts, x in Acts} GstarG[x,y] < ' char(39) fn_GstarG char(39) ';'])%--Read in matrix in AMPL
               % DEBUGGING: % %ampl{ni}.eval('printf {x in Acts, y in Acts}: "%.8g \n", GstarG[x,y] > ("GstarG_out_test.DAT"); ');
        % %--This matrix is too large for the Java heap space (unless it is increased). Write and read instead.

              
        ampl{ni}.eval('param RealGstarEab {q in Acts};');
        RealGstarEab = ampl{ni}.getParameter('RealGstarEab');
        RealGstarEab.setValues(cvar.RealGstarEab_wsum); %--Transfer values from Matlab to AMPL

        %--Regularization Matrix: Read in the diagonal, make the matrix of zeros, and fill in the diagonal.
        ampl{ni}.eval('param EyeGstarGdiag {q in Acts};');
        EyeGstarGdiag = ampl{ni}.getParameter('EyeGstarGdiag');
        EyeGstarGdiag.setValues(cvar.EyeGstarGdiag); %--Transfer values from Matlab to AMPL
        ampl{ni}.eval('param RegMat {x in Acts, y in Acts};');
        ampl{ni}.eval('let {x in Acts, y in Acts} RegMat[x,y] := 0;');
        ampl{ni}.eval('let {x in Acts, y in Acts: x==y} RegMat[x,y] := EyeGstarGdiag[x];');
        

        %--Solve
        ampl{ni}.eval('var costStep1 {x in Acts}  = sum {y in Acts} du[y]*(GstarG[x,y] + 10^(log10reg)*RegMat[x,y]);'); 
        ampl{ni}.eval('var costStep2 = sum {q in Acts}  (costStep1[q] + 2*RealGstarEab[q])*du[q];');  %--The cost function
        ampl{ni}.eval('minimize total_cost: costStep2; ')    
        %--The equations for reference:
        % ([u1; u2].' * (cvar.GstarG_wsum + 10^log10reg*maxDiagGstarG*eye(size(cvar.GstarG_wsum))) + cvar.RealGstarEab_wsum.') * [u1; u2] <= maxContrast
        % dDMvec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;
        
        
        % %  THIS WAY DOESN'T WORK!!!
%         ampl{ni}.eval('var costStep1 {q in Acts};'); 
%         ampl{ni}.eval('var costStep2;');  %--The cost function
%         ampl{ni}.eval('minimize total_cost: costStep2; ')     
%         ampl{ni}.eval('subject to costStep1_def {x in Acts}:  costStep1[x] = sum {y in Acts} du[y]*(GstarG[x,y] + 10^(log10reg)*RegMat[x,y]);');
%         ampl{ni}.eval('subject to costStep2_def:  costStep2 = sum {q in Acts}  (costStep1[q] + 2*RealGstarEab[q])*du[q];')

        ampl{ni}.eval('option solver gurobi;');
        %ampl{ni}.eval('option gurobi_options "outlev=1 threads=24";');
        %ampl{ni}.eval('option gurobi_options "outlev=1 presolve=0 lpmethod=2 barhomogeneous=1 crossover=0 threads=4";');
        ampl{ni}.eval(['option gurobi_options "outlev=1 presolve=0 lpmethod=2 barhomogeneous=1 crossover=0 threads=' num2str(floor(26/size(vals_list,2))) '";']);

        
        ampl{ni}.eval('solve;');
        

        %--Retrieve results
        totalcost = ampl{ni}.getObjective('total_cost');
        fprintf('Objective is: %f\n', totalcost.value);

        % Get the values of the variable Buy in a dataframe object
        du_temp1 = ampl{ni}.getVariable('du');
        du_temp2 = du_temp1.getValues;
        du_temp3 = du_temp2.val; whos a
        du_out = cell2mat(du_temp3); whos b
        % Print them
        % du_out

        %--Close AMPL
        ampl{ni}.reset();
        
        ampl{ni} = 0; %--try to overwrite to get rid of it for real
        
        
        %--Assign the output commands
%         if(any(mp.dm_ind==1)); dDM1V = zeros(mp.dm1.Nact,mp.dm1.Nact);  dDM1V(mp.dm1.act_ele) = du_out(u_dm_guide==1);  dDM1V_store(:,:,ni) = dDM1V; end
%         if(any(mp.dm_ind==2)); dDM2V = zeros(mp.dm2.Nact,mp.dm2.Nact);  dDM2V(mp.dm2.act_ele) = du_out(u_dm_guide==2);  dDM2V_store(:,:,ni) = dDM2V; end
%         if(any(mp.dm_ind==8)); dDM8V = du_out(u_dm_guide==8);  dDM8V_store(u_dm_guide==8,ni) = dDM8V; end
%         if(any(mp.dm_ind==9)); dDM9V = du_out(u_dm_guide==9);  dDM9V_store(u_dm_guide==9,ni) = dDM9V; end
        if(any(mp.dm_ind==1)); dDM.dDM1V = zeros(mp.dm1.Nact,mp.dm1.Nact);  dDM.dDM1V(mp.dm1.act_ele) = du_out(u_dm_guide==1);   end
        if(any(mp.dm_ind==2)); dDM.dDM2V = zeros(mp.dm2.Nact,mp.dm2.Nact);  dDM.dDM2V(mp.dm2.act_ele) = du_out(u_dm_guide==2);   end
        if(any(mp.dm_ind==8)); dDM.dDM8V = zeros(mp.dm8.NactTotal,1);  dDM.dDM8V(mp.dm8.act_ele) = du_out(u_dm_guide==8);  end
        if(any(mp.dm_ind==9)); dDM.dDM9V = zeros(mp.dm9.NactTotal,1);  dDM.dDM9V(mp.dm9.act_ele) = du_out(u_dm_guide==9);  end        
        
        %--Update the DM commands by adding the new control signal
        if(any(mp.dm_ind==1))
            mp.dm1.dV = dDM.dDM1V;
            mp.dm1.V = mp.dm1.V + mp.dm1.dV; 
        end
        if(any(mp.dm_ind==2))
            mp.dm2.dV = dDM.dDM2V;
            mp.dm2.V = mp.dm2.V + mp.dm2.dV; 
        end
        if(any(mp.dm_ind==8))
            mp.dm8.dV = dDM.dDM8V(:);
            mp.dm8.V = mp.dm8.V + mp.dm8.dV(:);
        end
        if(any(mp.dm_ind==9))
            mp.dm9.dV = dDM.dDM9V;
            mp.dm9.V(:) = mp.dm9.V(:) + mp.dm9.dV(:);
        end
        
        %%--Take images to empirically check contrast at that setting
        mp.flagParfor = false; %--Don't generate the image with parfor within parfor
        Itotal = falco_get_summed_image(mp);
        InormAvg = mean(Itotal(mp.F4.corr.maskBool));
        %Inorm_list(ni) = mean(Itotal(mp.F4.corr.maskBool));
        
    end %--END OF NESTED FUNCTION