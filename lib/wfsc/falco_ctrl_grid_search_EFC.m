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
% - Modified on 2019-09-26 by A.J. Riggs to handle DM actuator constraints
% outside this function in a more user-robust way.
% - Modified on 2019-06-25 by A.J. Riggs to pass out tied actuator pairs. 
% - Modified on 2018-07-24 to use Erkin's latest controller strategy.
% - Modified on 2018-02-06 by A.J. Riggs to be parallelized with parfor.
%   Required calling a new function. 
% - Modified by A.J. Riggs on October 11, 2017 to allow easier mixing of
%   which DMs are used and to also do a grid search over the gain of the 
%   overall DM command. 
% - Modified from hcil_ctrl_checkMuEmp.m by A.J. Riggs on August 31, 2016
% - Created at Princeton on 19 Feb 2015 by A.J. Riggs

function [dDM,cvarOut] = falco_ctrl_grid_search_EFC(mp,cvar)

    %--STEPS:
    % Step 0: [Done at begging of WFSC loop function] For this iteration, remove un-used DMs from the controller by changing mp.dm_ind value. 
    % Step 1: If re-linearizing this iteration, empirically find the best regularization value.
    % Step 2: For this iteration in the schedule, replace the imaginary part of the regularization with the latest "optimal" regularization
    % Step 3: Compute the EFC command to use.
    
    %% Initializations    
    vals_list = allcomb(mp.ctrl.log10regVec,mp.ctrl.dmfacVec).'; %--dimensions: [2 x length(mp.ctrl.muVec)*length(mp.ctrl.dmfacVec) ]
    Nvals = max(size(vals_list,2));
    Inorm_list = zeros(Nvals,1);

    % Temporarily store computed DM commands so that the best one does not have to be re-computed
    if(any(mp.dm_ind==1)); dDM1V_store = zeros(mp.dm1.Nact,mp.dm1.Nact,Nvals); end
    if(any(mp.dm_ind==2)); dDM2V_store = zeros(mp.dm2.Nact,mp.dm2.Nact,Nvals); end
    if(any(mp.dm_ind==5)); dDM5V_store = zeros(mp.dm5.Nact,mp.dm5.Nact,Nvals); end
    if(any(mp.dm_ind==8)); dDM8V_store = zeros(mp.dm8.NactTotal,Nvals); end
    if(any(mp.dm_ind==9)); dDM9V_store = zeros(mp.dm9.NactTotal,Nvals); end

    %% Empirically find the regularization value giving the best contrast
    
    %--Loop over all the settings to check empirically
    if(mp.flagParfor) %--Parallelized
        parfor ni = 1:Nvals
            [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar);
            %--delta voltage commands
            if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
            if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
            if(any(mp.dm_ind==5)); dDM5V_store(:,:,ni) = dDM_temp.dDM5V; end
            if(any(mp.dm_ind==8)); dDM8V_store(:,ni) = dDM_temp.dDM8V; end
            if(any(mp.dm_ind==9)); dDM9V_store(:,ni) = dDM_temp.dDM9V; end
        end
    else %--Not Parallelized
        for ni = Nvals:-1:1
            [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar);
            %--delta voltage commands
            if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
            if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
            if(any(mp.dm_ind==5)); dDM5V_store(:,:,ni) = dDM_temp.dDM5V; end
            if(any(mp.dm_ind==8)); dDM8V_store(:,ni) = dDM_temp.dDM8V; end
            if(any(mp.dm_ind==9)); dDM9V_store(:,ni) = dDM_temp.dDM9V; end
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
    [cvarOut.cMin,indBest] = min(Inorm_list(:));
    %--delta voltage commands
    if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest); end
    if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest); end
    if(any(mp.dm_ind==5)); dDM.dDM5V = dDM5V_store(:,:,indBest); end
    if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V_store(:,indBest); end
    if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V_store(:,indBest); end

    cvarOut.log10regUsed = vals_list(1,indBest);
    dmfacBest = vals_list(2,indBest);
    if(mp.ctrl.flagUseModel)
        fprintf('Model-based grid search gives log10reg, = %.1f,\t dmfac = %.2f\t   gives %4.2e contrast.\n',cvarOut.log10regUsed, dmfacBest, cvarOut.cMin)
    else
        fprintf('Empirical grid search gives log10reg, = %.1f,\t dmfac = %.2f\t   gives %4.2e contrast.\n',cvarOut.log10regUsed, dmfacBest, cvarOut.cMin)
    end
    
    %% Plot the grid search results
    if(mp.flagPlot)
        if(length(mp.ctrl.dmfacVec)==1)
            figure(499); semilogy(mp.ctrl.log10regVec,Inorm_list,'-bd','Linewidth',3)
            title('Line Search EFC','Fontsize',20,'Interpreter','Latex');
            xlabel('log10(regularization)','Fontsize',20,'Interpreter','Latex');
            ylabel('log10(Inorm)','Fontsize',20,'Interpreter','Latex');
            set(gca,'Fontsize',20); set(gcf,'Color',[1 1 1]); grid on;
            drawnow;
        elseif(length(mp.ctrl.dmfacVec)>1)
            figure(499); imagesc(mp.ctrl.log10regVec,mp.ctrl.dmfacVec,reshape(log10(Inorm_list),[length(mp.ctrl.dmfacVec),length(mp.ctrl.log10regVec)])); 
            ch = colorbar; axis xy tight;
            title('Grid Search EFC','Fontsize',20,'Interpreter','Latex');
            xlabel('log10(regularization)','Fontsize',20,'Interpreter','Latex');
            ylabel('Proportional Gain','Fontsize',20,'Interpreter','Latex');
            ylabel(ch,'log10(Inorm)','Fontsize',20,'Interpreter','Latex');
            set(gca,'Fontsize',20); set(gcf,'Color',[1 1 1]);
            drawnow;
        end
    end

end %--END OF FUNCTION