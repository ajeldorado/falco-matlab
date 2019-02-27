% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to
%
% -The algorithm in this code is mostly the same as electric field conjugation 
% (EFC) as described by Give'on et al. SPIE 2011. However, in this
% algorithm, it is the total DM stroke being minimized instead of the delta
% stroke in addition to the intensity.
%
%
% INPUTS:
% cvar.
%      dmfac = proportional gain factor
%      gamma = weight on total DM stroke (for actual DMs, # 1 and 2)
%      maxDiagGstarG
%
%
%
% REVISION HISTORY:
% -Created on 2018-11-06 by A.J. Riggs.

function [dDM,cvar,InormAvg] = falco_ctrl_total_stroke_minimization(mp,cvar)

%%
cvar.dampFac = mp.ctrl.dampFac;
cvar.dmfac = mp.ctrl.dmfac;
cvar.gamma = mp.ctrl.gamma;
cvar.gamma9 = mp.ctrl.gamma9; %--weight on total DM stroke for DM9.

    %% Make the vector of input, total control commands
    if(any(mp.dm_ind==1));  u1 = mp.dm1.V(mp.dm1.act_ele);  else;  u1 = [];  end
    if(any(mp.dm_ind==2));  u2 = mp.dm2.V(mp.dm2.act_ele);  else;  u2 = [];  end
    if(any(mp.dm_ind==5));  u5 = mp.dm5.V(mp.dm5.act_ele);  else;  u5 = [];  end
    if(any(mp.dm_ind==8));  u8 = mp.dm8.V(mp.dm8.act_ele);  else;  u8 = [];  end
    if(any(mp.dm_ind==9));  u9 = mp.dm9.V(mp.dm9.act_ele);  else;  u9 = [];  end
    u = [u1; u2; u5; u8; u9]; %--column vector
    %Nele = length(u);
    
    %--Get the indices of each DM's command within the full command
    if(any(mp.dm_ind==1));  u1dummy = 1*ones(mp.dm1.Nele,1);  else;  u1dummy = [];  end
    if(any(mp.dm_ind==2));  u2dummy = 2*ones(mp.dm2.Nele,1);  else;  u2dummy = [];  end
    if(any(mp.dm_ind==5));  u5dummy = 5*ones(mp.dm5.Nele,1);  else;  u5dummy = [];  end
    if(any(mp.dm_ind==8));  u8dummy = 8*ones(mp.dm8.Nele,1);  else;  u8dummy = [];  end
    if(any(mp.dm_ind==9));  u9dummy = 9*ones(mp.dm9.Nele,1);  else;  u9dummy = [];  end
    u_dm_guide = [u1dummy; u2dummy; u5dummy; u8dummy; u9dummy];

    %%
    
    GstarGdiag = diag(cvar.GstarG_wsum);
    cvar.maxDiagGstarG = max(GstarGdiag(u_dm_guide==1 | u_dm_guide==2));

    
    %%
%--Build up the diagonal of the regularization matrix. Allow different
%regularization values for different DMs
EyeDiag = [];
u_weight = [];
for idm=1:numel(mp.dm_ind)
    dm_index = mp.dm_ind(idm);
    if(dm_index==1)
        dm_reg = cvar.gamma;
        Nele = mp.dm1.Nele;
    elseif(dm_index==2)
        dm_reg = cvar.gamma;
        Nele = mp.dm2.Nele;
    elseif(dm_index==8)
        dm_reg = cvar.gamma8;
        Nele = mp.dm8.Nele;
    elseif(dm_index==9)
        dm_reg = cvar.gamma9;
        Nele = mp.dm9.Nele;
    else
        error('falco_ctrl_total_stroke_minimization.m: Controller weight not defined for DM %d.',dm_index);
    end
    
    EyeDiag = [EyeDiag; (dm_reg*cvar.maxDiagGstarG)*ones(Nele,1)];
    u_weight = [u_weight; (dm_reg*cvar.maxDiagGstarG)*ones(Nele,1)];
end


%% Least-squares solution. Different from EFC because of term u_weight.*u

dDMvec = -cvar.dmfac*(diag(EyeDiag) + cvar.GstarG_wsum)\(cvar.RealGstarEab_wsum + cvar.dampFac*u_weight.*u);
% dDMvec = -cvar.dmfac*(cvar.gamma*diag(EyeDiag) + cvar.GstarG_wsum)\(cvar.RealGstarEab_wsum + u_weight.*u);
    

%%
%--Initialize delta DM commands
if(any(mp.dm_ind==1)); dDM1V = zeros(mp.dm1.Nact); end
if(any(mp.dm_ind==2)); dDM2V = zeros(mp.dm2.Nact); end
if(any(mp.dm_ind==5)); dDM5V = zeros(mp.dm5.Nact); end
if(any(mp.dm_ind==8)); dDM8V = zeros(mp.dm8.NactTotal,1); end
if(any(mp.dm_ind==9)); dDM9V = zeros(mp.dm9.NactTotal,1); end

%% Parse the command vector to split the parts for DMs 1 & 2.
startIndex = 0; % Initialize starting index of the command vector

%--DM1
if(any(mp.dm_ind==1)) %--DM1
    dDM1V(mp.dm1.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm1.Nele); % Parse the command vector to get component for DM1
    dDM1V = dDM1V*mp.dm1.weight; %--Re-scale correctly based on the DM's weighting
    dDM1Vmax = max(abs(dDM1V(:))); % Store max absolute deviation value for later
    startIndex = startIndex + mp.dm1.Nele; % Set for next DM
end

%--DM2
if(any(mp.dm_ind==2)) %--DM2
    dDM2V(mp.dm2.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm2.Nele);
    dDM2V = dDM2V*mp.dm2.weight; %--Re-scale correctly based on the DM's weighting
    dDM2Vmax = max(abs(dDM2V(:)));
    startIndex = startIndex + mp.dm2.Nele; % Set for next DM
end

%--DM5
if(any(mp.dm_ind==5)) %--DM5
    dDM5V(mp.dm5.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm5.Nele);
    dDM5V = dDM5V*mp.dm_weights(5); %--Re-scale correctly based on the DM's weighting
    dDM5Vmax = max(abs(dDM5V(:)));
    startIndex = startIndex + mp.dm5.Nele; % Set for next DM
end

%--DM8
if(any(mp.dm_ind==8))  %--DM8
    dDM8V(mp.dm8.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm8.Nele);
    dDM8V = dDM8V*mp.dm8.weight; %--Re-scale correctly based on the DM's weighting
    startIndex = startIndex + mp.dm8.Nele; % Set for next DM
end

%--DM9
if(any(mp.dm_ind==9))  %--DM9
    dDM9V(mp.dm9.act_ele) = dDMvec(startIndex+1:startIndex+mp.dm9.Nele);
    dDM9V = dDM9V*mp.dm9.weight; %--Re-scale correctly based on the DM's weighting
    startIndex = startIndex + mp.dm9.Nele; % Set for next DM
end

 
%% Add the delta command to the command

%--Save starting point for each delta command to be added to.
if(any(mp.dm_ind==1)); DM1Vnom = mp.dm1.V; end
if(any(mp.dm_ind==2)); DM2Vnom = mp.dm2.V; end
if(any(mp.dm_ind==5)); DM5Vnom = mp.dm5.V; end
if(any(mp.dm_ind==8)); DM8Vnom = mp.dm8.V(:); end
if(any(mp.dm_ind==9)); DM9Vnom = mp.dm9.V; end

if(any(mp.dm_ind==1)); mp.dm1.V = DM1Vnom + dDM1V; end
if(any(mp.dm_ind==2)); mp.dm2.V = DM2Vnom + dDM2V; end
if(any(mp.dm_ind==5)); mp.dm5.V = DM5Vnom + dDM5V;  end
if(any(mp.dm_ind==8)); mp.dm8.V(:) = DM8Vnom(:) + dDM8V(:);  end
if(any(mp.dm_ind==8)); mp.dm9.V = DM9Vnom + dDM9V;  end

if(any(mp.dm_ind==1)); dDM.dDM1V = dDM1V; end
if(any(mp.dm_ind==2)); dDM.dDM2V = dDM2V; end
if(any(mp.dm_ind==5)); dDM.dDM5V = dDM5V; end
if(any(mp.dm_ind==8)); dDM.dDM8V = dDM8V(:); end
if(any(mp.dm_ind==9)); dDM.dDM9V = dDM9V(:); end


%% Take images to empirically check contrast at that mu value
Itotal = falco_get_summed_image(mp);
InormAvg = mean(Itotal(mp.Fend.corr.maskBool));
cvar.cMin = InormAvg;

cvar.log10regUsed = -10; %--Dummy value



end
