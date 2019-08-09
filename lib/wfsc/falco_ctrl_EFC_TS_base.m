% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
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
%
% REVISION HISTORY: 
% - Modified on 2019-07-08 by A.J. Riggs from falco_ctrl_EFC_base.m to
% falco_ctrl_EFC_total_stroke_base.m to have term for total stroke.
% - Modified on 2019-04-23 by A.J. Riggs to have an option for a
%   model-based grid search.
% - Modified on 2019-02-13 by A.J. Riggs to use falco_ctrl_setup.m and
%   falco_ctrl_wrapup.m to save a bunch of space.
% - Modified on 2018-11-12 by A.J. Riggs to clean up code, especially to 
%   remove the large commented-out blocks.
% - Modified on 2018-07-24 by A.J. Riggs to switch from the Lagrange multiplier to the Tikhonov regularization.
% - Modified on 2018-02-06 by A.J. Riggs to be parallelized with parfor.
%   Called by a higher function. 
% - Modified by A.J. Riggs on October 11, 2017 to allow easier mixing of
%   which DMs are used and to also do a grid search over the gain of the 
%   overall DM command. 
% - Modified from hcil_ctrl_checkMuEmp.m by A.J. Riggs on August 31, 2016
% - Created at Princeton on 19 Feb 2015 by A.J. Riggs

%--Return values:
%  Measured average normalized intensity
%  DM commands

function [InormAvg,dDM] = falco_ctrl_EFC_TS_base(ni,vals_list,mp,cvar)


%% Initializations
log10reg = vals_list(1,ni); %--Lagrange multiplier
dmfac = vals_list(2,ni); %--Scaling factor for entire DM command

%--Save starting point for each delta command to be added to.
%--Get the indices of each DM's command vector within the single concatenated command vector
cvar = falco_ctrl_setup(mp,cvar);

%--Static weights on total stroke
gammaVec = mp.ctrl.gamma*ones(size(diag(cvar.EyeGstarGdiag)));
if(any(mp.dm_ind==8)); gammaVec(cvar.uLegend==8) = mp.ctrl.gamma8; end
if(any(mp.dm_ind==9)); gammaVec(cvar.uLegend==9) = mp.ctrl.gamma9; end


%% Least-squares solution with regularization:
% duVec = -dmfac*( diag(10^(log10reg) + gammaVec)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\(cvar.RealGstarEab_wsum + cvar.EyeNorm*gammaVec.*cvar.uVec);
duVec = -dmfac*( 10^(log10reg)*diag(cvar.EyeGstarGdiag) + diag(gammaVec.*cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\(cvar.RealGstarEab_wsum + cvar.EyeNorm*gammaVec.*cvar.uVec);

%% Parse the command vector by DM and assign the output commands
[mp,dDM] = falco_ctrl_wrapup(mp,cvar,duVec);

%% Take images and compute average intensity in dark hole

if(mp.ctrl.flagUseModel) %--Perform a model-based grid search using the compact model
    Itotal = falco_get_expected_summed_image(mp,cvar);
    InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
else %--Perform an empirical grid search with actual images from the testbed or full model
    Itotal = falco_get_summed_image(mp);
    InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
end
        
end %--END OF FUNCTION