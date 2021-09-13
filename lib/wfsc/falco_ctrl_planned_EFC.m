% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Peform EFC according to a pre-defined (planned) scheme.
%
% INPUTS
% ------
% mp : structure of model parameters
% cvar : structure of controller variables
%
% OUTPUTS
% -------
% dDM : structure of the delta control commands separated by DM number.
%       Also contains the updated array of tied actuator pairs
% cvar : structure of controller variables

function [dDM, cvar] = falco_ctrl_planned_EFC(mp, cvar)

    %--STEPS:
    % Step 0: [Done at begging of WFSC loop function] For this iteration, remove un-used DMs from the controller by changing mp.dm_ind value. 
    % Step 1: If re-linearizing this iteration, empirically find the best regularization value.
    % Step 2: For this iteration in the schedule, replace the imaginary part of the regularization with the latest "optimal" regularization
    % Step 3: Compute the EFC command to use.
    
    %% Initializations    
    vals_list = allcomb(mp.ctrl.log10regVec,mp.ctrl.dmfacVec).'; %--dimensions: [2 x length(mp.ctrl.muVec)*length(mp.ctrl.dmfacVec) ]
    Nvals = max(size(vals_list,2));
    Inorm_list = zeros(Nvals,1);

    %--Use these to temporarily store computed DM commands so that the best one does not have to be re-computed
    if(any(mp.dm_ind==1)); dDM1V_store = zeros(mp.dm1.Nact,mp.dm1.Nact,Nvals); end
    if(any(mp.dm_ind==2)); dDM2V_store = zeros(mp.dm2.Nact,mp.dm2.Nact,Nvals); end
    if(any(mp.dm_ind==8)); dDM8V_store = zeros(mp.dm8.NactTotal,Nvals); end
    if(any(mp.dm_ind==9)); dDM9V_store = zeros(mp.dm9.NactTotal,Nvals); end

    %% Step 1: Empirically find the "optimal" regularization value (if told to for this iteration)
    
    if(any(mp.gridSearchItrVec==cvar.Itr))
        
        %--Loop over all the settings to check empirically
        ImCube = zeros(mp.Fend.Neta, mp.Fend.Nxi, Nvals);
        if mp.flagParfor && (mp.flagSim || mp.ctrl.flagUseModel) %--Parallelized
            parfor ni = 1:Nvals
                [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar);
                %--delta voltage commands
                if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
                if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
                if(any(mp.dm_ind==8)); dDM8V_store(:,ni) = dDM_temp.dDM8V; end
                if(any(mp.dm_ind==9)); dDM9V_store(:,ni) = dDM_temp.dDM9V; end
                ImCube(:, :, ni) = dDM_temp.Itotal;
            end
        else %--Not Parallelized
            for ni = Nvals:-1:1
                [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar);
                %--delta voltage commands
                if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
                if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
                if(any(mp.dm_ind==8)); dDM8V_store(:,ni) = dDM_temp.dDM8V; end
                if(any(mp.dm_ind==9)); dDM9V_store(:,ni) = dDM_temp.dDM9V; end
                ImCube(:, :, ni) = dDM_temp.Itotal;
            end
        end

        %--Print out results to the command line
        fprintf('Scaling factor:\t')
        for ni=1:Nvals;  fprintf('%.2f\t\t', vals_list(2,ni) );  end

        fprintf('\nlog10reg:\t');
        for ni=1:Nvals;  fprintf('%.1f\t\t',vals_list(1,ni));  end

        fprintf('\nInorm:  \t')
        for ni=1:Nvals;  fprintf('%.2e\t',Inorm_list(ni));  end
        fprintf('\n')

        %--Find the best scaling factor and Lagrange multiplier pair based on the best contrast.
        [cvar.cMin,indBest] = min(Inorm_list(:));
        cvar.latestBestlog10reg = vals_list(1,indBest);
        cvar.latestBestDMfac = vals_list(2,indBest);
        cvar.Im = ImCube(:, :, indBest);
        if(mp.ctrl.flagUseModel)
            fprintf('Model-based grid search expects log10reg, = %.1f,\t dmfac = %.2f\t   gives %4.2e normalized intensity.\n',cvar.latestBestlog10reg, cvar.latestBestDMfac, cvar.cMin)
        else
            fprintf('Empirical grid search finds log10reg, = %.1f,\t dmfac = %.2f\t   gives %4.2e normalized intensity.\n',cvar.latestBestlog10reg, cvar.latestBestDMfac, cvar.cMin)
        end
    end
    
    %% Skip steps 2 and 3 if the schedule for this iteration is just to use the "optimal" regularization AND if the grid search was performed this iteration
    
    if( (imag(mp.ctrl.log10regSchedIn(cvar.Itr)) ~= 0) && (real(mp.ctrl.log10regSchedIn(cvar.Itr)) == 0) && any(mp.gridSearchItrVec==cvar.Itr) )
        %--delta voltage commands
        if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest); end
        if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest); end
        if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V_store(:,indBest); end
        if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V_store(:,indBest); end
        
        log10regSchedOut = cvar.latestBestlog10reg;
    else
    
        %% Step 2: For this iteration in the schedule, replace the imaginary part of the regularization with the latest "optimal" regularization
        if( imag(mp.ctrl.log10regSchedIn(cvar.Itr)) ~= 0)
            log10regSchedOut = cvar.latestBestlog10reg + real(mp.ctrl.log10regSchedIn(cvar.Itr));
        else
            log10regSchedOut = mp.ctrl.log10regSchedIn(cvar.Itr);
        end

        %% Step 3: Compute the EFC command to use.
        ni = 1;
        if(isfield(cvar,'latestBestDMfac')==false)
            cvar.latestBestDMfac = 1;
        end
        vals_list = [log10regSchedOut; cvar.latestBestDMfac];
        
        [cvar.cMin,dDM] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar);
        cvar.Im = dDM.Itotal;
        dDM = rmfield(dDM, 'Itotal'); % reduce amount of memory used since image moved to cvar structure
        if(mp.ctrl.flagUseModel)
            fprintf('Model says scheduled log10reg = %.1f\t gives %4.2e contrast.\n',log10regSchedOut,cvar.cMin)
        else
            fprintf('Scheduled log10reg = %.1f\t gives %4.2e contrast.\n',log10regSchedOut,cvar.cMin)
        end
    end
    
    cvar.log10regUsed = log10regSchedOut;
       
end %--END OF FUNCTION
