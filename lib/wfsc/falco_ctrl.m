% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Wrapper function for the different controller functions.
%
% INPUTS
% ------
% mp : structure of model parameters
% cvar : structure of controller variables
% jacStruct : structure containing the Jacobians
%
% OUTPUTS
% -------
% mp : structure of model parameters
% cvar : structure of controller variables


function [mp, cvar] = falco_ctrl(mp, cvar, jacStruct)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Control Algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    fprintf('Using the Jacobian to make other matrices...'); tic;

    %--Compute matrices for linear control with regular EFC
    cvar.GstarG_wsum  = zeros(cvar.NeleAll, cvar.NeleAll); 
    cvar.RealGstarEab_wsum = zeros(cvar.NeleAll, 1);
    Eest = cvar.Eest;
    
    jacStruct = falco_apply_spatial_weighting_to_Jacobian(mp, jacStruct);
    
    for iMode = 1:mp.jac.Nmode

        Gstack = [jacStruct.G1(:,:,iMode), jacStruct.G2(:,:,iMode), jacStruct.G3(:,:,iMode), jacStruct.G4(:,:,iMode), jacStruct.G5(:,:,iMode), jacStruct.G6(:,:,iMode), jacStruct.G7(:,:,iMode), jacStruct.G8(:,:,iMode), jacStruct.G9(:,:,iMode)];

        %--Square matrix part stays the same if no re-linearization has occurrred. 
        cvar.GstarG_wsum  = cvar.GstarG_wsum  + mp.jac.weights(iMode) * real(Gstack'*Gstack); 

        %--The G^*E part changes each iteration because the E-field changes.
        if mp.jac.minimizeNI
            modvar.whichSource = 'star';
            modvar.sbpIndex = mp.jac.sbp_inds(iMode);
            modvar.zernIndex = mp.jac.zern_inds(iMode);
            modvar.starIndex = mp.jac.star_inds(iMode);
            Eunocculted = model_compact(mp, modvar, 'nofpm');
            [~, indPeak] = max(abs(Eunocculted(:)));
            Epeak = Eunocculted(indPeak);
            Eest(:, iMode) = cvar.Eest(:, iMode) / Epeak;           
        end
        
        iStar = mp.jac.star_inds(iMode);
        Eweighted = mp.WspatialVec(:, iStar) .* Eest(:, iMode); %--Apply 2-D spatial weighting to E-field in dark hole pixels.
        cvar.RealGstarEab_wsum = cvar.RealGstarEab_wsum + mp.jac.weights(iMode)*real(Gstack'*Eweighted); %--Apply the Jacobian weights and add to the total.

    end
    clear Gstack Eweighted % save RAM
    
    %--Make the regularization matrix. (Define only the diagonal here to save RAM.)
    cvar.EyeGstarGdiag = max(diag(cvar.GstarG_wsum ))*ones(cvar.NeleAll,1);
    cvar.EyeNorm = max(diag(cvar.GstarG_wsum ));
    fprintf(' done. Time: %.3f\n',toc);

    %--Call the Controller Function
    fprintf('Control beginning ...\n'); tic
    switch lower(mp.controller)

        %--Established, conventional controllers
        case{'plannedefc'} %--EFC regularization is scheduled ahead of time
            [dDM,cvar] = falco_ctrl_planned_EFC(mp,cvar);

        case{'gridsearchefc'}  %--Empirical grid search of EFC. Scaling factor for DM commands too.
            [dDM,cvar] = falco_ctrl_grid_search_EFC(mp,cvar);
            
        
        %--Experimental controllers
        case{'plannedefcts'} %--EFC regularization is scheduled ahead of time. total stroke also minimized
            [dDM,cvar] = falco_ctrl_planned_EFC_TS(mp,cvar);
            
        case{'plannedefccon'} %--Constrained-EFC regularization is scheduled ahead of time
            [dDM,cvar] = falco_ctrl_planned_EFCcon(mp,cvar);
            
        case{'sm-cvx'} %--Constrained & bounded stroke minimization using CVX. The quadratic cost function is solved directly CVX.
            cvar.dummy = 1;
            [dDM,cvar] = falco_ctrl_SM_CVX(mp,cvar);
            
        case{'tsm'}
            cvar.dummy = 1;
            [dDM,cvar] = falco_ctrl_total_stroke_minimization(mp,cvar); 
            
    end
    fprintf(' done. Time: %.3f sec\n',toc);

    mp = falco_ctrl_update_dm_commands(mp, dDM);
    
end
