% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_compact(mp, DM, modvar)
%--Blind model used by the estimator and controller
%  Does not include unknown aberrations/errors that are in the full model.
%
% REVISION HISTORY:
% --------------
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
% -modvar = structure of model variables
%
%
% OUTPUTS:
% -Eout
%  -> Note: When computing the control Jacobian, Eout is a structure
%  containing the Jacobians. Otherwise, it is just the E-field at the
%  second focal plane.
%
% modvar structure fields (4):
% -sbpIndex
% -wpsbpIndex
% -whichSource
% -flagGenMat


function Eout = model_compact(mp, DM, modvar)

%--Save a lot of RAM
if(any(DM.dm_ind==1)); DM.dm1.compact = rmfield(DM.dm1.compact,'inf_datacube'); end
if(any(DM.dm_ind==2)); DM.dm2.compact = rmfield(DM.dm2.compact,'inf_datacube'); end

%--Select the type of coronagraph
switch mp.coro 
        
    case{'LC','DMLC','APLC'} %--optional apodizer, FPM with/without phase contribution, and LS.
        Eout = model_compact_LC(mp, DM, modvar);  
       
    case{'SPLC','FLC'} %--Optional apodizer, binary-amplitude FPM with outer diaphragm, LS
        Eout = model_compact_SPLC(mp, DM, modvar);
            
    case{'vortex','Vortex','VC','AVC'} %--Optional apodizer, vortex FPM, LS
        Eout = model_compact_VC(mp, DM, modvar);      
        
%     case{'SPC','APP','APC'} %--Pupil-plane mask only
%         Eout = model_compact_APC(mp, DM, modvar);             

    otherwise
        disp('ERROR: CASE NOT RECOGNIZED IN model_compact.m');        
end
    

end % End of function


    
