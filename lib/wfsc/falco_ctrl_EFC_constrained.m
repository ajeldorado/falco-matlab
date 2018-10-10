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
% -Modified on 2018-02-06 by A.J. Riggs to be parallelized with parfor.
%   Required calling a new function. 
% -Modified by A.J. Riggs on October 11, 2017 to allow easier mixing of
%   which DMs are used and to also do a grid search over the gain of the 
%   overall DM command. 
% -Modified from hcil_ctrl_checkMuEmp.m by A.J. Riggs on August 31, 2016
% -Created at Princeton on 19 Feb 2015 by A.J. Riggs

function [dDM,cvar] = falco_ctrl_EFC_constrained(cp, cvar, mp)

    Itr = cvar.Itr ;

    if(any(mp.dm_ind==9))
        vals_list = allcomb(mp.ctrl.muVec,mp.ctrl.dm9regfacVec).'; %--dimensions: [2 x length(mp.ctrl.muVec)*length(mp.ctrl.dm9regfacVec) ]
    else
        vals_list = mp.ctrl.muVec; %--dimensions: [1 x length(mp.ctrl.muVec) ]
    end
    Nvals = max(size(vals_list,2));
    Inorm_list = zeros(Nvals,1);    
    
    
%     if(any(mp.dm_ind==9))
%         vals_list = allcomb(mp.ctrl.muVec,mp.ctrl.dmfacVec,mp.ctrl.dm9regfacVec).'; %--dimensions: [3 x length(mp.ctrl.muVec)*length(mp.ctrl.dmfacVec)*length(mp.ctrl.dm9regfacVec) ]
%     else
%         vals_list = allcomb(mp.ctrl.muVec,mp.ctrl.dmfacVec).'; %--dimensions: [2 x length(mp.ctrl.muVec)*length(mp.ctrl.dmfacVec) ]
%     end
%     Nvals = max(size(vals_list,2));
%     Inorm_list = zeros(Nvals,1);


    % Temporarily store computed DM commands so that the best one does not have to be re-computed
    if(any(mp.dm_ind==1)); dDM1V_store = zeros(mp.dm1.Nact,mp.dm1.Nact,Nvals); end
    if(any(mp.dm_ind==2)); dDM2V_store = zeros(mp.dm2.Nact,mp.dm2.Nact,Nvals); end
    if(any(mp.dm_ind==8)); dDM8V_store = zeros(mp.dm8.NactTotal,Nvals); end
    if(any(mp.dm_ind==9)); dDM9V_store = zeros(mp.dm9.NactTotal,Nvals); end


%     if(mp.flagParfor) %--Parallelized %--DOES NOT WORK WITH CVX. THE CVX PROGRAMS OVERWRITE EACH OTHER.
%         parfor ni = 1:Nvals
%             [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_constrained_func(ni,vals_list,DM, cvar);%(ni,vals_list,DM, cp, cvar, mp);
%             if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
%             if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
%             if(any(mp.dm_ind==9)); dDM9V_store(:,ni) = dDM_temp.dDM9V; end
%         end
%     else %--Not Parallelized
        for ni = 1:Nvals
            [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_constrained_func(ni,vals_list, cvar, mp);
            if(any(mp.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
            if(any(mp.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
            if(any(mp.dm_ind==8)); dDM8V_store(:,ni) = dDM_temp.dDM8V; end
            if(any(mp.dm_ind==9)); dDM9V_store(:,ni) = dDM_temp.dDM9V; end
        end
%     end


    %--Print out results to the command line if more than one used
if(numel(mp.ctrl.muVec)>1)

    fprintf('\nMu:\t\t');
    for ni=1:Nvals
        fprintf('%.2e\t',vals_list(1,ni))
    end

    fprintf('\nInorm:  \t')
    for ni=1:Nvals
        fprintf('%.2e\t',Inorm_list(ni))
    end
    fprintf('\n')

end



    %--Find the best Lagrange multiplier value based on the best contrast.
    [cvar.cMin,indBest] = min(Inorm_list(:));

    if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest); end
    if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest); end
    if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V_store(:,indBest); end
    if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V_store(:,indBest); end
  
    muBest = vals_list(1,indBest);
    fprintf('Empirically chosen mu = %.2e\t   gives %4.2e contrast.\n',muBest,cvar.cMin)    

    %cp.muBest(Itr) = vals_list(1,indBest);
    %fprintf('Empirically chosen mu = %.2e\t   gives %4.2e contrast.\n',cp.muBest(Itr),cvar.cMin)    


end %--END OF FUNCTION
