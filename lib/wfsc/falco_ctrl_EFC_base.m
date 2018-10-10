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
% -Modified on 2018-07-24 by A.J. Riggs to switch from the Lagrange multiplier to the Tikhonov regularization.
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
function [InormAvg,dDM] = falco_ctrl_EFC_base(ni,vals_list,mp,cvar)

% Itr = cvar.Itr ;

log10reg = vals_list(1,ni); %--Lagrange multiplier
dmfac = vals_list(2,ni); %--Scaling factor for entire DM command

%--Save starting point for each delta command to be added to.
if(any(mp.dm_ind==1)); DM1Vnom = mp.dm1.V; end
if(any(mp.dm_ind==2)); DM2Vnom = mp.dm2.V; end
if(any(mp.dm_ind==8)); DM8Vnom = mp.dm8.V(:); end
if(any(mp.dm_ind==9)); DM9Vnom = mp.dm9.V; end

% % %--Diagonal of the Weighted Regularization Matrix
% % EyeGstarGdiag = [];
% % maxDiagGstarG = max(diag(cvar.GstarG_wsum));
% % for idm=1:numel(mp.dm_ind)
% %     dm_index = mp.dm_ind(idm);
% %     dm_weight = 1; %mp.dm_weights(dm_index);
% %     if(any(mp.dm_ind==9)); dm_weight = dm9regfac*dm_weight; end
% %     EyeGstarGdiag = [EyeGstarGdiag; maxDiagGstarG*dm_weight*ones(cvar.NeleVec(idm),1)];
% % end

%--Least-squares solution:
dDMvec = -dmfac*(10^(log10reg)*diag(cvar.EyeGstarGdiag) + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;
% dDMvec = -dmfac*(diag(EyeGstarGdiag)/mu + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;
% % dDMvec = -dmfac*(diag(cvar.EyeGstarGdiag)/mu + cvar.GstarG_wsum)\cvar.RealGstarEab_wsum;

%--Initialize delta DM commands
if(any(mp.dm_ind==1)); dDM1V = zeros(mp.dm1.Nact); end
if(any(mp.dm_ind==2)); dDM2V = zeros(mp.dm2.Nact); end
if(any(mp.dm_ind==8)); dDM8V = zeros(mp.dm8.NactTotal,1); end
if(any(mp.dm_ind==9)); dDM9V = zeros(mp.dm9.NactTotal,1); end

%% Parse the command vector to split the parts for DMs 1 & 2.
startIndex = 0; % Initialize starting index of the command vector

%--DM1
if(any(mp.dm_ind==1)) %--DM1
    dDM1V(mp.dm1.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm1.Nele); % Parse the command vector to get component for DM1
    dDM1V = dDM1V*mp.dm_weights(1); %--Re-scale correctly based on the DM's weighting
    dDM1Vmax = max(abs(dDM1V(:))); % Store max absolute deviation value for later
    startIndex = startIndex + mp.dm1.Nele; % Set for next DM
end

%--DM2
if(any(mp.dm_ind==2)) %--DM2
    dDM2V(mp.dm2.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm2.Nele);
    dDM2V = dDM2V*mp.dm_weights(2); %--Re-scale correctly based on the DM's weighting
    dDM2Vmax = max(abs(dDM2V(:)));
    startIndex = startIndex + mp.dm2.Nele; % Set for next DM
end

%--DM8
if(any(mp.dm_ind==8))  %--DM8
    dDM8V(mp.dm8.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm8.Nele);
    dDM8V = dDM8V*mp.dm_weights(8); %--Re-scale correctly based on the DM's weighting
    startIndex = startIndex + mp.dm8.Nele; % Set for next DM

    %--Rescale the voltage range if it goes too high (and becomes highly nonlinear)
    if( max(abs(dDM8V(:))) > mp.dm8.dVmax )  
        dDM8V = dDM8V*(mp.dm8.dVmax/max(abs(dDM8V(:)))); 
    end
    
    %--Apply voltage range restrictions
    DM8Vtemp = DM8Vnom(:) + dDM8V(:);
    DM8Vtemp(DM8Vtemp < mp.dm8.Vmin) = mp.dm8.Vmin; 
    DM8Vtemp(DM8Vtemp > mp.dm8.Vmax) = mp.dm8.Vmax; 
%     DM8Vtemp(find(DM8Vtemp < mp.dm8.Vmin)) = mp.dm8.Vmin; 
%     DM8Vtemp(find(DM8Vtemp > mp.dm8.Vmax)) = mp.dm8.Vmax; 
    
    dDM8V(:) = DM8Vtemp(:) - DM8Vnom(:);
    mp.dm8.V(:) = DM8Vnom(:) + dDM8V(:);
    
end

%--DM9
if(any(mp.dm_ind==9))  %--DM9
    dDM9V(mp.dm9.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm9.Nele);
    dDM9V = dDM9V*mp.dm_weights(9); %--Re-scale correctly based on the DM's weighting
    startIndex = startIndex + mp.dm9.Nele; % Set for next DM
    
    %--Rescale the voltage range if it goes too high (and becomes highly nonlinear)
    if( max(abs(dDM9V(:))) > mp.maxAbsdV )  
        dDM9V = dDM9V*(mp.maxAbsdV/max(abs(dDM9V(:)))); 
    end
    
    %--Apply voltage range restrictions
    DM9Vtemp = DM9Vnom + dDM9V;
    DM9Vtemp(find(DM9Vtemp < mp.dm9.Vmin)) = mp.dm9.Vmin; 
    DM9Vtemp(find(DM9Vtemp > mp.dm9.Vmax)) = mp.dm9.Vmax; 
    
    dDM9V = DM9Vtemp - DM9Vnom;
    mp.dm9.V = DM9Vnom + dDM9V;
    
end

%% Rescale the DM1 and DM2 voltage ranges equally if either goes too high (and becomes highly nonlinear)  
if(any(mp.dm_ind==1) && any(mp.dm_ind==2)) %--DMs 1 and 2 used
    if( (dDM1Vmax > mp.maxAbsdV) || (dDM2Vmax > mp.maxAbsdV) )
        dDM1V = dDM1V*(mp.maxAbsdV/max([dDM1Vmax,dDM2Vmax]));   
        dDM2V = dDM2V*(mp.maxAbsdV/max([dDM1Vmax,dDM2Vmax]));   
    end
elseif(any(mp.dm_ind==1)) %--DM1,but no DM2
    if(dDM1Vmax > mp.maxAbsdV)
        dDM1V = dDM1V*(mp.maxAbsdV/dDM1Vmax); 
    end
elseif(any(mp.dm_ind==2)) %--DM2, but no DM1
    if(dDM2Vmax > mp.maxAbsdV)
        dDM2V = dDM2V*(mp.maxAbsdV/dDM2Vmax); 
    end
end
    
if(any(mp.dm_ind==1)) 
    %--Apply voltage range restrictions
    DM1Vtemp = DM1Vnom + dDM1V;
    DM1Vtemp(find(DM1Vtemp < -mp.dm1.maxAbsV)) = -mp.dm1.maxAbsV ; 
    DM1Vtemp(find(DM1Vtemp >  mp.dm1.maxAbsV)) =  mp.dm1.maxAbsV ; 
    
    dDM1V = DM1Vtemp - DM1Vnom;
    mp.dm1.V = DM1Vnom + dDM1V;
end

if(any(mp.dm_ind==2)) 
    %--Apply voltage range restrictions
    DM2Vtemp = DM2Vnom + dDM2V;
    DM2Vtemp(find(DM2Vtemp < -mp.dm2.maxAbsV)) = -mp.dm2.maxAbsV ; 
    DM2Vtemp(find(DM2Vtemp >  mp.dm2.maxAbsV)) =  mp.dm2.maxAbsV ; 
    
    dDM2V = DM2Vtemp - DM2Vnom;
    mp.dm2.V = DM2Vnom + dDM2V; 
end


%% Take images to empirically check contrast at that mu value
Itotal = falco_get_summed_image(mp);
InormAvg = mean(Itotal(mp.F4.corr.maskBool));
        

if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V; end
if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V; end
if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V(:); end
if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V(:); end

end %--END OF FUNCTION
