% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function jac = model_Jacobian(mp, DM, modvar)
%--Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%
% REVISION HISTORY:
% --------------
% Modified on 2018-03-06 by A.J. Riggs to be more compatible with parfor.
% Modified on 2017-11-13 by A.J. Riggs to be compatible with parfor. Added
% this function as an extra wrapper layer.
% Modified on 2017-11-09 by A.J. Riggs from model_compact.m to
%   model_Jacobian.m
% Modified on 2017-10-17 by A.J. Riggs to have model_compact.m be a wrapper. All the 
%  actual compact models have been moved to sub-routines for clarity.
% Modified on 19 June 2017 by A.J. Riggs to use lower resolution than the
%   full model.
% model_compact.m - 18 August 2016: Modified from hcil_model.m
% hcil_model.m - 18 Feb 2015: Modified from HCIL_model_lab_BB_v3.m
% ---------------
%
% INPUTS:
% -mp = structure of model parameters
% -DM = structure of DM settings
% -vals_list = structure containing combinations of:
%     -tsi = index of the pair of sub-bandpass index and tip/tilt offset index
%     -whichDM = which DM number
%
% OUTPUTS:
% -Jac = Jacobian for the specified DM and specified T/T-wavelength pair

function Jac = model_Jacobian_middle_layer(mp, DM, vals_list,ii)

    tsi = vals_list(1,ii); %--index for tip-tilt-&-sub-bandpass-pair
    whichDM = vals_list(2,ii); %--number of the specified DM

    switch mp.coro 
        
        case{'LC','DMLC','APLC'} %--pupil+DMs, occulting spot FPM, and LS.
            Jac = model_Jacobian_LC(mp, DM, tsi, whichDM); 
            
        case{'SPLC','FLC'} %--Optional apodizer, binary-amplitude FPM with outer diaphragm, LS
            Jac  = model_Jacobian_SPLC(mp, DM, tsi, whichDM); 
            
        case{'vortex','Vortex','VC','AVC'} %--Optional apodizer, vortex FPM, LS
            Jac  = model_Jacobian_VC(mp, DM, tsi, whichDM); 
            
        %case{'SPC','APP','APC'} %--Pupil-plane apodizer is only coronagraphic mask
            %Jac  = model_Jacobian_APC(mp, DM, tsi, whichDM); 
            
        otherwise
            error('model_Jacobian_middle_layer: CASE NOT RECOGNIZED IN model_Jacobian.m');        
    end    

end
