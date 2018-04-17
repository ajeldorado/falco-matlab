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

function [dDM,cp,cvar] = falco_ctrl_EFC(DM, cp, cvar, mp)

    Itr = cvar.Itr ;


    vals_list = allcomb(cp.muVec,cp.dmfacVec).'; %--dimensions: [2 x length(cp.muVec)*length(cp.dmfacVec) ]
    
    Nvals = max(size(vals_list,2));
    Inorm_list = zeros(Nvals,1);


    % Temporarily store computed DM commands so that the best one does not have to be re-computed
    if(any(DM.dm_ind==1)); dDM1V_store = zeros(DM.dm1.Nact,DM.dm1.Nact,Nvals); end
    if(any(DM.dm_ind==2)); dDM2V_store = zeros(DM.dm2.Nact,DM.dm2.Nact,Nvals); end


    if(mp.flagParfor) %--Parallelized
        parfor ni = 1:Nvals
            [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_func(ni,vals_list,DM, cp, cvar, mp);
            if(any(DM.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
            if(any(DM.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
        end
    else %--Not Parallelized
        for ni = 1:Nvals
            [Inorm_list(ni),dDM_temp] = falco_ctrl_EFC_func(ni,vals_list,DM, cp, cvar, mp);
            if(any(DM.dm_ind==1)); dDM1V_store(:,:,ni) = dDM_temp.dDM1V; end
            if(any(DM.dm_ind==2)); dDM2V_store(:,:,ni) = dDM_temp.dDM2V; end
        end
    end


    %--Print out results to the command line
    fprintf('Scaling factor:\t')
    for ni=1:Nvals
        fprintf('%.2f\t\t', vals_list(2,ni) )
    end

    fprintf('\nMu:\t\t');
    for ni=1:Nvals
        fprintf('%.2e\t',vals_list(1,ni))
    end

    fprintf('\nInorm:  \t')
    for ni=1:Nvals
        fprintf('%.2e\t',Inorm_list(ni))
    end
    fprintf('\n')
    


    %--Find the best scaling factor and Lagrange multiplier pair based on the best contrast.
    [cvar.cMin,indBest] = min(Inorm_list(:));

    if(any(DM.dm_ind==1)); dDM.dDM1V = dDM1V_store(:,:,indBest); end
    if(any(DM.dm_ind==2)); dDM.dDM2V = dDM2V_store(:,:,indBest); end



    cp.muBest(Itr) = vals_list(1,indBest);
    cp.dmfacBest(Itr) = vals_list(2,indBest);
    fprintf('Empirically chosen mu = %.2e,\t dmfac = %.2f\t   gives %4.2e contrast.\n',cp.muBest(Itr),cp.dmfacBest(Itr),cvar.cMin)
    

end %--END OF FUNCTION
