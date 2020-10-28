% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function for regularized linear least-squares control (EFC).
% -This function computes the DM commands and new Inorm for one set of these parameters:
%  a) a scalar coefficient for the regularization matrix
%  b) a scalar gain for the final DM command.
%
% -This code is based on electric field conjugation (EFC) as described 
% by Give'on et al. SPIE 2011.
%
%--Return values:
%  Measured average normalized intensity
%  DM commands

function [InormAvg,dDM] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar)


%% Initializations
log10reg = vals_list(1,ni); %--Lagrange multiplier
dmfac = vals_list(2,ni); %--Scaling factor for entire DM command

%--Save starting point for each delta command to be added to.
%--Get the indices of each DM's command vector within the single concatenated command vector
cvar = falco_ctrl_setup(mp,cvar);

%% Least-squares solution with regularization:
duVec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;

%% Parse the command vector by DM and assign the output commands
[mp,dDM] = falco_ctrl_wrapup(mp,cvar,duVec);

%% Take images and compute average intensity in dark hole

if(mp.ctrl.flagUseModel) %--Perform a model-based grid search using the compact model
    if(mp.flagFiber == false)
        Itotal = falco_get_expected_summed_image(mp,cvar);
        InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
        dDM.Itotal = Itotal;
    end
else %--Perform an empirical grid search with actual images from the testbed or full model
    if(mp.flagFiber == false)
        Itotal = falco_get_summed_image(mp);
        InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
        dDM.Itotal = Itotal;
    else
        IfiberTotal = falco_get_summed_image_fiber(mp);
        if(mp.flagLenslet)
            InormAvg = mean(max(max(IfiberTotal)));
        else
            InormAvg = mean(IfiberTotal(mp.Fend.corr.maskBool));
        end

    end
end
        
end %--END OF FUNCTION