% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function for regularized linear least-squares control (EFC).
% -This function performs an empirical grid search over these parameters:
%  a) a scalar coefficient for the regularization matrix
%  b) a scalar gain for the final DM command.
%
% -This code is based on electric field conjugation (EFC) as described 
% by Give'on et al. SPIE 2011.
%
%
% REVISION HISTORY:
% -Modified on 2018-07-24 to use Erkin's latest controller strategy.
% -Modified on 2018-02-06 by A.J. Riggs to be parallelized with parfor.
%   Required calling a new function. 
% -Modified by A.J. Riggs on October 11, 2017 to allow easier mixing of
%   which DMs are used and to also do a grid search over the gain of the 
%   overall DM command. 
% -Modified from hcil_ctrl_checkMuEmp.m by A.J. Riggs on August 31, 2016
% -Created at Princeton on 19 Feb 2015 by A.J. Riggs

function [dDM,cvar] = falco_ctrl_planned_EFC(mp, cvar)

    %--STEPS:
    % Step 0: [Done at begging of WFSC loop function] For this iteration, remove un-used DMs from the controller by changing mp.dm_ind value. 
    % Step 1: If re-linearizing this iteration, empirically find the best regularization value.
    % Step 2: For this iteration in the schedule, replace the imaginary part of the regularization with the latest "optimal" regularization
    % Step 3: Compute the EFC command to use.
    
    %% Initializations    
    if cvar.Itr<mp.aux.ItrDump 
        vals_list = allcomb(linspace(0,3,7),mp.ctrl.dmfacVec).';
    else
        vals_list = allcomb(mp.ctrl.log10regVec,mp.ctrl.dmfacVec).'; %--dimensions: [2 x length(mp.ctrl.muVec)*length(mp.ctrl.dmfacVec) ]
    end
    Nvals = max(size(vals_list,2));
    if mp.aux.flagOmega==1 && cvar.Itr>=mp.aux.firstOmegaItr
        NvalsOmega = 5;
        valsOmega_list = [linspace(mp.aux.omegaMin,mp.aux.omegaMax,NvalsOmega)];
    else
        valsOmega_list = [-inf];
    end
    NvalsOmega = numel(valsOmega_list);
    Inorm_list = zeros(Nvals,NvalsOmega);
    thput_list = zeros(Nvals,NvalsOmega);

    %--Use these to temporarily store computed DM commands so that the best one does not have to be re-computed
    if(any(mp.dm_ind==1)); dDM1V_store = zeros(mp.dm1.Nact,mp.dm1.Nact,Nvals,NvalsOmega); end
    if(any(mp.dm_ind==2)); dDM2V_store = zeros(mp.dm2.Nact,mp.dm2.Nact,Nvals,NvalsOmega); end
    if(any(mp.dm_ind==5)); dDM5V_store = zeros(mp.dm5.Nact,mp.dm5.Nact,Nvals,NvalsOmega); end
    if(any(mp.dm_ind==8)); dDM8V_store = zeros(mp.dm8.NactTotal,Nvals,NvalsOmega); end
    if(any(mp.dm_ind==9)); dDM9V_store = zeros(mp.dm9.NactTotal,Nvals,NvalsOmega); end

    %% Step 1: Empirically find the "optimal" regularization value (if told to for this iteration)
    
    if(any(mp.gridSearchItrVec==cvar.Itr))
        
        %--Loop over all the settings to check empirically
        if(mp.flagParfor) %--Parallelized
            parfor ni = 1:Nvals
                for nj = 1:NvalsOmega
                [Inorm_list(ni,nj),thput_list(ni,nj),dDM_temp] = falco_ctrl_EFC_base(ni,vals_list,nj,valsOmega_list,mp,cvar);
                if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni,nj) = dDM_temp.dDM1V; end
                if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni,nj) = dDM_temp.dDM2V; end
                if(any(mp.dm_ind==5)); dDM5V_store(:,:,ni,nj) = dDM_temp.dDM5V; end
                if(any(mp.dm_ind==8)); dDM8V_store(:,ni,nj) = dDM_temp.dDM8V; end
                if(any(mp.dm_ind==9)); dDM9V_store(:,ni,nj) = dDM_temp.dDM9V; end
                end
             end
        else %--Not Parallelized
            for ni = 1:Nvals
                for nj = 1:NvalsOmega
                [Inorm_list(ni,nj),thput_list(ni,nj),dDM_temp] = falco_ctrl_EFC_base(ni,vals_list,nj,valsOmega_list,mp,cvar);
                if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni,nj) = dDM_temp.dDM1V; end
                if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni,nj) = dDM_temp.dDM2V; end
                if(any(mp.dm_ind==5)); dDM5V_store(:,:,ni,nj) = dDM_temp.dDM5V; end
                if(any(mp.dm_ind==8)); dDM8V_store(:,ni,nj) = dDM_temp.dDM8V; end
                if(any(mp.dm_ind==9)); dDM9V_store(:,ni,nj) = dDM_temp.dDM9V; end
                end
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
%         [cvar.cMin,indBest] = min(Inorm_list(:));
        [cvar.cMin,indBest] = min(Inorm_list(:)./(thput_list(:).^2));
        
        [indBest,indBestOmega,indBestRegDM9] = ind2sub(size(Inorm_list),indBest);
        mp.aux.omega = valsOmega_list(indBestOmega);
        
        val = vals_list(1,indBest)-mp.aux.betaMinusOne;
        if ismember(val,vals_list)
            indBest=find(vals_list(1,:)==val);
            cvar.latestBestlog10reg = vals_list(1,indBest);
            cvar.omegaUsed = mp.aux.omega;
            cvarOut.regDM9Used = vals_list_dm9(indBestRegDM9);
            cvar.cMin = Inorm_list(indBest);
            dmfacBest = vals_list(2,indBest);
            cvar.latestBestDMfac = vals_list(2,indBest);
            if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest,indBestOmega); end
            if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest,indBestOmega); end
            if(any(mp.dm_ind==5)); dDM.dDM5V = dDM5V_store(:,:,indBest,indBestOmega); end
            if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V_store(:,indBest,indBestOmega); end
            if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V_store(:,indBest,indBestOmega); end
        else
            vals_listaux = [val;vals_list(2,indBest)];
            [Inorm,thput,dDM_temp] = falco_ctrl_EFC_base(1,vals_listaux,mp,cvar);
            cvar.latestBestlog10reg = val;
            cvar.cMin = Inorm;
            dmfacBest = vals_list(2,indBest);
            cvar.latestBestDMfac = vals_list(2,indBest);
            if(any(mp.dm_ind==1)); dDM.dDM1V = dDM_temp.dDM1V; end
            if(any(mp.dm_ind==2)); dDM.dDM2V = dDM_temp.dDM2V; end
            if(any(mp.dm_ind==5)); dDM.dDM5V = dDM_temp.dDM5V; end
            if(any(mp.dm_ind==8)); dDM.dDM8V = dDM_temp.dDM8V; end
            if(any(mp.dm_ind==9)); dDM.dDM9V = dDM_temp.dDM9V; end
        end

%         cvar.latestBestlog10reg = vals_list(1,indBest);
%         cvar.latestBestDMfac = vals_list(2,indBest);
        
        %cp.muBest(cvar.Itr) = vals_list(1,indBest);
        %cp.dmfacBest(cvar.Itr) = vals_list(2,indBest);
        fprintf('Empirical grid search gives log10reg, = %.1f,\t dmfac = %.2f\t   gives %4.2e contrast.\n',cvar.latestBestlog10reg,cvar.latestBestDMfac,cvar.cMin)

    end
    
    %% Skip steps 2 and 3 if the schedule for this iteration is just to use the "optimal" regularization AND if the grid search was performed this iteration
    
    if( (imag(mp.ctrl.log10regSchedIn(cvar.Itr)) ~= 0) && (real(mp.ctrl.log10regSchedIn(cvar.Itr)) == 0) && any(mp.gridSearchItrVec==cvar.Itr) )
    
        if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest,indBestOmega); end
        if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest,indBestOmega); end
        if(any(mp.dm_ind==5)); dDM.dDM5V = dDM5V_store(:,:,indBest,indBestOmega); end
        if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V_store(:,indBest,indBestOmega); end
        if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V_store(:,indBest,indBestOmega); end
        
        log10regSchedOut = cvar.latestBestlog10reg;
    else
    
        %% Step 2: For this iteration in the schedule, replace the imaginary part of the regularization with the latest "optimal" regularization
        if( imag(mp.ctrl.log10regSchedIn(cvar.Itr)) ~= 0)
            log10regSchedOut = cvar.latestBestlog10reg + real(mp.ctrl.log10regSchedIn(cvar.Itr));
        else
            log10regSchedOut = mp.ctrl.log10regSchedIn(cvar.Itr);
        end
    %     fprintf('log10reg = %.1f \n',log10regSchedOut )
    %     fprintf('mp.dm_ind = ['); for jj = 1:length(mp.dm_ind); fprintf(' %d',mp.dm_ind(jj)); end; fprintf(' ]\n');


        %% Step 3: Compute the EFC command to use.
        ni = 1;
        if(isfield(cvar,'latestBestDMfac')==false)
            cvar.latestBestDMfac = 1;
        end
        vals_list = [log10regSchedOut; cvar.latestBestDMfac];
<<<<<<< Updated upstream
        
        [cvar.cMin,dDM] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar);
        fprintf('Scheduled log10reg = %.1f\t gives %4.2e contrast.\n',log10regSchedOut,cvar.cMin)
%         fprintf('Scheduled values of: \tlog10reg = %.1f,\t dmfac = %.2f\t   gives %4.2e contrast.\n',log10regSchedOut,cvar.latestBestDMfac,cvar.cMin)

=======
        [cvar.cMin,thput,dDM] = falco_ctrl_EFC_base(ni,vals_list,indBestOmega,valsOmega_list,indBestRegDM9,vals_list_dm9,mp,cvar);
%         [cvar.cMin,dDM] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar);
        if(mp.ctrl.flagUseModel)
            fprintf('Model says scheduled log10reg = %.1f\t gives %4.2e contrast.\n',log10regSchedOut,cvar.cMin)
        else
            fprintf('Scheduled log10reg = %.1f\t gives %4.2e contrast.\n',log10regSchedOut,cvar.cMin)
        end
>>>>>>> Stashed changes
    end
    
    cvar.log10regUsed = log10regSchedOut;
    
    
end %--END OF FUNCTION
