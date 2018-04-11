% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
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
% -Modified on 2018-02-06 by A.J. Riggs to be parallelized with parfor.
%   Called by a higher function. 
% -Modified by A.J. Riggs on October 11, 2017 to allow easier mixing of
%   which DMs are used and to also do a grid search over the gain of the 
%   overall DM command. 
% -Modified from hcil_ctrl_checkMuEmp.m by A.J. Riggs on August 31, 2016
% -Created at Princeton on 19 Feb 2015 by A.J. Riggs


%--Return values:
%  Measured average normalized intensity
%  DM commands
function [InormAvg,dDM] = falco_ctrl_EFC_func(ni,vals_list,DM, cp, cvar, mp)

% Itr = cvar.Itr ;

mu = vals_list(1,ni); %--Lagrange multiplier
dmfac = vals_list(2,ni); %--Scaling factor for entire DM command



%--Save starting point for each delta command to be added to.
if(any(DM.dm_ind==1)); DM1Vnom = DM.dm1.V; end
if(any(DM.dm_ind==2)); DM2Vnom = DM.dm2.V; end


%--Diagonal of the Weighted Regularization Matrix
EyeGstarGdiag = [];
maxDiagGstarG = max(diag(cvar.GstarG_wsum));
for idm=1:numel(DM.dm_ind)
    dm_index = DM.dm_ind(idm);
    dm_weight = DM.dm_weights(dm_index);
    EyeGstarGdiag = [EyeGstarGdiag; maxDiagGstarG*dm_weight*ones(cvar.NeleVec(idm),1)];
end

%--Least-squares solution:
dDMvec = -dmfac*(diag(EyeGstarGdiag)/mu + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;
% dDMvec = -dmfac*(diag(cvar.EyeGstarGdiag)/mu + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;

%--Initialize delta DM commands
if(any(DM.dm_ind==1)); dDM1V = zeros(DM.dm1.Nact); end
if(any(DM.dm_ind==2)); dDM2V = zeros(DM.dm2.Nact); end

%--Parse the command vector to split the parts for each DM.
startIndex = 0; % Initialize starting index of the command vector
if(any(DM.dm_ind==1)) %--DM1
    dDM1V(DM.dm1.act_ele) = dDMvec(startIndex+1:startIndex+DM.dm1.Nele);
    dDM1Vmax = max(abs(dDM1V(:)));
    startIndex = startIndex + DM.dm1.Nele; % Set for next DM
end

if(any(DM.dm_ind==2)) %--DM2
    dDM2V(DM.dm2.act_ele) = dDMvec(startIndex+1:startIndex+DM.dm2.Nele);
    dDM2Vmax = max(abs(dDM2V(:)));
    startIndex = startIndex + DM.dm2.Nele; % Set for next DM
end


%--Rescale the DM1 and DM2 voltage ranges equally if either goes too high (and becomes highly nonlinear)  
if(any(DM.dm_ind==1) && any(DM.dm_ind==2)) %--DMs 1 and 2 used
    if( (dDM1Vmax > DM.maxAbsdV) || (dDM2Vmax > DM.maxAbsdV) );
        dDM1V = dDM1V*(DM.maxAbsdV/max([dDM1Vmax,dDM2Vmax]));   
        dDM2V = dDM2V*(DM.maxAbsdV/max([dDM1Vmax,dDM2Vmax]));   
    end
elseif(any(DM.dm_ind==1)) %--DM1,but no DM2
    if(dDM1Vmax > DM.maxAbsdV);
        dDM1V = dDM1V*(DM.maxAbsdV/dDM1Vmax); 
    end
elseif(any(DM.dm_ind==2)) %--DM2, but no DM1
    if(dDM2Vmax > DM.maxAbsdV);
        dDM2V = dDM2V*(DM.maxAbsdV/dDM2Vmax); 
    end
end
    

if(any(DM.dm_ind==1)); 
    %--Apply voltage range restrictions
    DM1Vtemp = DM1Vnom + dDM1V;
    DM1Vtemp(find(DM1Vtemp < -DM.dm1.maxAbsV)) = -DM.dm1.maxAbsV ; 
    DM1Vtemp(find(DM1Vtemp >  DM.dm1.maxAbsV)) =  DM.dm1.maxAbsV ; 
    
    dDM1V = DM1Vtemp - DM1Vnom;
    DM.dm1.V = DM1Vnom + dDM1V;
end
if(any(DM.dm_ind==2)); 
    %--Apply voltage range restrictions
    DM2Vtemp = DM2Vnom + dDM2V;
    DM2Vtemp(find(DM2Vtemp < -DM.dm2.maxAbsV)) = -DM.dm2.maxAbsV ; 
    DM2Vtemp(find(DM2Vtemp >  DM.dm2.maxAbsV)) =  DM.dm2.maxAbsV ; 
    
    dDM2V = DM2Vtemp - DM2Vnom;
    DM.dm2.V = DM2Vnom + dDM2V; 
end


% Take images to empirically check contrast at that mu value
Itotal = falco_get_summed_image(mp, DM);
InormAvg = mean(Itotal(mp.F4.full.corr.inds));
        

if(any(DM.dm_ind==1)); dDM.dDM1V = dDM1V; end
if(any(DM.dm_ind==2)); dDM.dDM2V = dDM2V; end

end %--END OF FUNCTION
