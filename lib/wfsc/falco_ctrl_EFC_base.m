% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute EFC commands and take images for the chosen regularization and
% scale factor values.
%
% INPUTS
% ------
% iList : index to use in valsList
% valsList : list of all parameter value combinations. made by allcomb
% mp : structure of model parameters
% cvar : structure of controller variables
%
% OUTPUTS
% -------
% InormMean : 2-D intensity image averaged over all subbands
% dDM : structure of the delta control commands separated by DM number.

function [InormMean, dDM] = falco_ctrl_EFC_base(iList, valsList, mp, cvar)

log10reg = valsList(1, iList);
scaleFactor = valsList(2, iList);

cvar = falco_ctrl_setup(mp, cvar);

duVec = -scaleFactor*((10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum); % One-line EFC equation

[mp, dDM] = falco_ctrl_wrapup(mp, cvar, duVec);

%% Take images and compute average intensity in dark hole

if mp.ctrl.flagUseModel %--Perform model-based grid search using compact model
    if mp.flagFiber
        error('Fiber option not implemented for mp.ctrl.flagUseModel==true.')
    elseif strcmpi(mp.estimator, 'iefc')
        error("Cannot use mp.ctrl.flagUseModel=true when 'iefc' is mp.estimator.")
    else
        Itotal = falco_get_expected_summed_image(mp,cvar);
        InormMean = mean(Itotal(mp.Fend.corr.maskBool));
        dDM.Itotal = Itotal;
    end
    
else %--Perform an empirical grid search with actual images
    if ~mp.flagFiber
        Itotal = falco_get_summed_image(mp);
        InormMean = mean(Itotal(mp.Fend.corr.maskBool));
        dDM.Itotal = Itotal;
    else
%         IfiberTotal = falco_get_summed_image_fiber(mp);
        [~,IfiberTotal] = falco_get_summed_image(mp);
        dDM.Itotal = IfiberTotal;
        if(mp.flagLenslet)
            InormMean = mean(max(max(IfiberTotal)));
        else
            InormMean = IfiberTotal;
        end

    end
end
        
end %--END OF FUNCTION
