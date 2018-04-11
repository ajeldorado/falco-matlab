% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_full(mp, DM, modvar)
%--Full-knowledge optical model.
%    --> Not used by the estimator and controller.
%    --> Only used to create simulated intensity images.
%
% REVISION HISTORY:
% --------------
% Modified on 2017-10-17 by A.J. Riggs to have model_full.m be a wrapper. All the 
%  actual full models have been moved to sub-routines for clarity.
% model_full.m - Modified from hcil_simTestbed.m
% hcil_simTestbed.m - 18 Feb 2015: Modified from hcil_model.m. Includes
%  extra errors in the model to simulate the actual testbed for fake images.
%
% ---------------
% INPUTS:
% -mp = structure of model parameters
% -DM = structure of DM settings
% -modvar = structure of model variables
%
%
% OUTPUTS:
% -Eout
%
% modvar structure fields (4):
% -sbpIndex
% -wpsbpIndex
% -whichSource
% -flagGenMat


function Eout = model_full(mp, DM, modvar)

%--Save a lot of RAM (remove influence function data cubes, except for DM9
% which needs it.
if(any(DM.dm_ind==1)); DM.dm1 = rmfield(DM.dm1,'compact'); end
if(any(DM.dm_ind==2)); DM.dm2 = rmfield(DM.dm2,'compact'); end

%--Select the type of coronagraph
switch mp.coro 

    case{'LC','DMLC','APLC'} %--optional apodizer, occulting spot FPM, and LS.
        Eout = model_full_LC(mp, DM, modvar);
       
    case{'SPLC','FLC'} %--Optional apodizer, binary-amplitude FPM with outer diaphragm, LS
        Eout = model_full_SPLC(mp, DM, modvar);
            
    case{'vortex','Vortex','VC','AVC'} %--Optional apodizer, vortex FPM, LS
        Eout = model_full_VC(mp, DM, modvar);       
          
%     case{'SPC','APP','APC'} %--Pupil-plane mask only
%         Eout = model_full_APC(mp, DM, modvar);   
        
    otherwise
        fprintf('ERROR: CASE %s NOT RECOGNIZED IN models/model_full.m\n',mp.coro);        
                    
end

if(mp.useGPU)
    Eout = gather(Eout);
end

end % End of function


    
