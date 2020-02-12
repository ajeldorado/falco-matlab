% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function acting as a wrapper for the different controller functions.
%
% INPUTS:
% -mp = structure of model parameters
% -cvar = structure containing variables for the controller
% -jacStruct = structure containing the Jacobians
%
% OUTPUTS:
% -mp = structure of model parameters
% REVISION HISTORY:
% --------------
% Moved on 2019-09-30 by A.J. Riggs to be a separate file. Had tried
%   nesting in falco_wfsc_loop to save RAM but that did not work.
% Created by A.J. Riggs on 2018-10-04 by extracting material from falco_wfsc_loop.m.
% ---------------

function [mp,cvar] = falco_ctrl(mp,cvar,jacStruct)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Control Algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    fprintf('Using the Jacobian to make other matrices...'); tic;

    %--Compute matrices for linear control with regular EFC
    cvar.GstarG_wsum  = zeros(cvar.NeleAll,cvar.NeleAll); 
    cvar.RealGstarEab_wsum = zeros(cvar.NeleAll, 1);

    for im=1:mp.jac.Nmode

        Gstack = [jacStruct.G1(:,:,im), jacStruct.G2(:,:,im), jacStruct.G3(:,:,im), jacStruct.G4(:,:,im), jacStruct.G5(:,:,im), jacStruct.G6(:,:,im), jacStruct.G7(:,:,im), jacStruct.G8(:,:,im), jacStruct.G9(:,:,im)];

        %--Square matrix part stays the same if no re-linearization has occurrred. 
        cvar.GstarG_wsum  = cvar.GstarG_wsum  + mp.jac.weights(im)*real(Gstack'*Gstack); 

        %--The G^*E part changes each iteration because the E-field changes.
        Eweighted = mp.WspatialVec.*cvar.EfieldVec(:,im); %--Apply 2-D spatial weighting to E-field in dark hole pixels.
        cvar.RealGstarEab_wsum = cvar.RealGstarEab_wsum + mp.jac.weights(im)*real(Gstack'*Eweighted); %--Apply the Jacobian weights and add to the total.

    end
    clear GallCell Gstack Eweighted % save RAM
    
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
    
    %% Updates to DM commands

    %--Update the DM commands by adding the delta control signal
    if(any(mp.dm_ind==1));  mp.dm1.V = mp.dm1.V + dDM.dDM1V;  end
    if(any(mp.dm_ind==2));  mp.dm2.V = mp.dm2.V + dDM.dDM2V;  end
    if(any(mp.dm_ind==3));  mp.dm3.V = mp.dm3.V + dDM.dDM3V;  end
    if(any(mp.dm_ind==4));  mp.dm4.V = mp.dm4.V + dDM.dDM4V;  end
    if(any(mp.dm_ind==5));  mp.dm5.V = mp.dm5.V + dDM.dDM5V;  end
    if(any(mp.dm_ind==6));  mp.dm6.V = mp.dm6.V + dDM.dDM6V;  end
    if(any(mp.dm_ind==7));  mp.dm7.V = mp.dm7.V + dDM.dDM7V;  end
    if(any(mp.dm_ind==8));  mp.dm8.V = mp.dm8.V + dDM.dDM8V;  end
    if(any(mp.dm_ind==9));  mp.dm9.V = mp.dm9.V + dDM.dDM9V;  end

    %%--Save the delta from the previous command
    if(any(mp.dm_ind==1));  mp.dm1.dV = dDM.dDM1V;  end
    if(any(mp.dm_ind==2));  mp.dm2.dV = dDM.dDM2V;  end
    if(any(mp.dm_ind==3));  mp.dm3.dV = dDM.dDM3V;  end
    if(any(mp.dm_ind==4));  mp.dm4.dV = dDM.dDM4V;  end
    if(any(mp.dm_ind==5));  mp.dm5.dV = dDM.dDM5V;  end
    if(any(mp.dm_ind==6));  mp.dm6.dV = dDM.dDM6V;  end
    if(any(mp.dm_ind==7));  mp.dm7.dV = dDM.dDM7V;  end
    if(any(mp.dm_ind==8));  mp.dm8.dV = dDM.dDM8V;  end
    if(any(mp.dm_ind==9));  mp.dm9.dV = dDM.dDM9V;  end

end %--END OF FUNCTION falco_ctrl
