% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to
%
%
%
% REVISION HISTORY:
% -Created on 2018-11-06 by A.J. Riggs.

function [dDM,cvar,InormAvg] = falco_ctrl_SVD_truncation_EFC(mp, jacStruct,cvar)

keyboard

    %% Make the vector of input, total control commands
    if(any(mp.dm_ind==1));  u1 = mp.dm1.V(mp.dm1.act_ele);  else;  u1 = [];  end
    if(any(mp.dm_ind==2));  u2 = mp.dm2.V(mp.dm2.act_ele);  else;  u2 = [];  end
    if(any(mp.dm_ind==3));  u3 = mp.dm3.V(mp.dm3.act_ele);  else;  u3 = [];  end
    if(any(mp.dm_ind==4));  u4 = mp.dm4.V(mp.dm4.act_ele);  else;  u4 = [];  end
    if(any(mp.dm_ind==5));  u5 = mp.dm5.V(mp.dm5.act_ele);  else;  u5 = [];  end
    if(any(mp.dm_ind==6));  u6 = mp.dm6.V(mp.dm6.act_ele);  else;  u6 = [];  end
    if(any(mp.dm_ind==7));  u7 = mp.dm7.V(mp.dm7.act_ele);  else;  u7 = [];  end
    if(any(mp.dm_ind==8));  u8 = mp.dm8.V(mp.dm8.act_ele);  else;  u8 = [];  end
    if(any(mp.dm_ind==9));  u9 = mp.dm9.V(mp.dm9.act_ele);  else;  u9 = [];  end
    u = [u1; u2; u3; u4; u5; u6; u7; u8; u9]; %--column vector
    NeleAll = length(u);

    %--Get the indices of each DM's command within the full command
    if(any(mp.dm_ind==1));  u1dummy = 1*ones(mp.dm1.Nele,1);  else;  u1dummy = [];  end
    if(any(mp.dm_ind==2));  u2dummy = 2*ones(mp.dm2.Nele,1);  else;  u2dummy = [];  end
    if(any(mp.dm_ind==3));  u3dummy = 3*ones(mp.dm3.Nele,1);  else;  u3dummy = [];  end
    if(any(mp.dm_ind==4));  u4dummy = 4*ones(mp.dm4.Nele,1);  else;  u4dummy = [];  end
    if(any(mp.dm_ind==5));  u5dummy = 5*ones(mp.dm5.Nele,1);  else;  u5dummy = [];  end
    if(any(mp.dm_ind==6));  u6dummy = 6*ones(mp.dm6.Nele,1);  else;  u6dummy = [];  end
    if(any(mp.dm_ind==7));  u7dummy = 7*ones(mp.dm7.Nele,1);  else;  u7dummy = [];  end
    if(any(mp.dm_ind==8));  u8dummy = 8*ones(mp.dm8.Nele,1);  else;  u8dummy = [];  end
    if(any(mp.dm_ind==9));  u9dummy = 9*ones(mp.dm9.Nele,1);  else;  u9dummy = [];  end
    cvar.u_guide = [u1dummy; u2dummy; u3dummy;  u4dummy;  u5dummy; u6dummy;  u7dummy;  u8dummy; u9dummy];

    %--Save starting point for each delta command to be added to.
    if(any(mp.dm_ind==1)); cvar.DM1Vnom = mp.dm1.V; end
    if(any(mp.dm_ind==2)); cvar.DM2Vnom = mp.dm2.V; end
    if(any(mp.dm_ind==3)); cvar.DM3Vnom = mp.dm3.V; end
    if(any(mp.dm_ind==4)); cvar.DM4Vnom = mp.dm4.V; end
    if(any(mp.dm_ind==5)); cvar.DM5Vnom = mp.dm5.V; end
    if(any(mp.dm_ind==6)); cvar.DM6Vnom = mp.dm6.V; end
    if(any(mp.dm_ind==7)); cvar.DM7Vnom = mp.dm7.V; end
    if(any(mp.dm_ind==8)); cvar.DM8Vnom = mp.dm8.V(:); end
    if(any(mp.dm_ind==9)); cvar.DM9Vnom = mp.dm9.V; end

%%
    %--Make the Jacobian into one giant matrix
    Npix = size(cvar.EfieldVec,1);
    EvecAll = reshape(cvar.EfieldVec,[Npix*mp.Nsbp, 1]);
    Gall = zeros(Npix*mp.Nsbp,cvar.NeleAll);
    for im=1:mp.jac.Nmode
        Gall( (im-1)*Npix+1:im*Npix ,:) = [jacStruct.G1(:,:,im), jacStruct.G2(:,:,im), jacStruct.G3(:,:,im), jacStruct.G4(:,:,im), jacStruct.G5(:,:,im), jacStruct.G6(:,:,im), jacStruct.G7(:,:,im), jacStruct.G8(:,:,im), jacStruct.G9(:,:,im)];
    end
    
%     GstarG = Gall'*Gall;
%     alphaSq = max(diag(GstarG));
%     %%
%     EvecAll2 = zeros(Npix*mp.Nsbp,1);
%     for im=1:mp.jac.Nmode
%         EvecAll2((im-1)*Npix+1:im*Npix) = cvar.EfieldVec(:,im);
%     end
%     figure; plot(abs(EvecAll2-EvecAll));
%%
    %--Perform SVD
%     [U,S,V] = svd(Gall,0);
    [U,S,V] = svd(Gall);

    %%
    %--Perform truncation of singular values
%     mp.NItarget = 1e-10;
    
    Esvd = U'*EvecAll;
    Idh = sum(abs(EvecAll).^2);
    Ndh = size(Gall,1);
    
    
    IsvdSum = zeros(size(diag(S)));
    IdiffVec = zeros(size(diag(S)));
    for ii = 1:length(diag(S))
        IsvdSum(ii) = sum( abs(Esvd(1:ii)).^2 )/Ndh;
        IdiffVec(ii) = (Idh - sum( abs(Esvd(1:ii)).^2 ))/Ndh;
    end
    figure(21); semilogy(IsvdSum); set(gca,'Fontsize',20);
    figure(22); semilogy(IdiffVec); set(gca,'Fontsize',20);

    
    Idiff = 1;
    counter = 1; 
    while Idiff > mp.NItarget
        Idiff = (Idh - sum( abs(Esvd(1:counter)).^2 ))/Ndh;
        fprintf('%4d:   Idiff = %.3e\n',counter,Idiff);
        counter = counter + 1;
        
        %--Break if going on forever
        if(counter > length(diag(S)))
            disp('      *** No singular values truncated in falco_ctrl_SVD_truncated_EFC.m. ***       ')
            Idiff = mp.NItarget-1;
        end
    end
    counter = counter - 1;
    
    counter = 1000;
    
    %--Truncated singular values
    Strunc = 0*S;
    Sdiag = diag(S);
    Strunc(1:counter) = Sdiag(1:counter);
%     SdiagTrunc = Sdiag;
%     SdiagTrunc(counter+1:end) = 0;
    
    %--Control command
    duVec = real( V*Strunc'.^-1*Esvd );
%     dDMvec = V*diag(1/SdiagTrunc)*Esvd;


%% Parse the command vector by DM

%--Initialize delta DM commands
if(any(mp.dm_ind==1)); dDM.dDM1V = zeros(mp.dm1.Nact,mp.dm1.Nact); end
if(any(mp.dm_ind==2)); dDM.dDM2V = zeros(mp.dm2.Nact,mp.dm2.Nact); end
% if(any(mp.dm_ind==1)); dDM.dDM1V = zeros(mp.dm1.NactTotal,1); end
% if(any(mp.dm_ind==2)); dDM.dDM2V = zeros(mp.dm2.NactTotal,1); end
if(any(mp.dm_ind==3)); dDM.dDM3V = zeros(mp.dm3.NactTotal,1); end
if(any(mp.dm_ind==4)); dDM.dDM4V = zeros(mp.dm4.NactTotal,1); end
if(any(mp.dm_ind==5)); dDM.dDM5V = zeros(mp.dm5.NactTotal,1); end
if(any(mp.dm_ind==6)); dDM.dDM6V = zeros(mp.dm6.NactTotal,1); end
if(any(mp.dm_ind==7)); dDM.dDM7V = zeros(mp.dm7.NactTotal,1); end
if(any(mp.dm_ind==8)); dDM.dDM8V = zeros(mp.dm8.NactTotal,1); end
if(any(mp.dm_ind==9)); dDM.dDM9V = zeros(mp.dm9.NactTotal,1); end

%--Parse the command vector by DM
if(any(mp.dm_ind==1));  dDM.dDM1V(mp.dm1.act_ele) = mp.dm_weights(1)*duVec(cvar.u_guide==1);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==2));  dDM.dDM2V(mp.dm2.act_ele) = mp.dm_weights(2)*duVec(cvar.u_guide==2);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==3));  dDM.dDM3V(mp.dm3.act_ele) = mp.dm_weights(3)*duVec(cvar.u_guide==3);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==4));  dDM.dDM4V(mp.dm4.act_ele) = mp.dm_weights(4)*duVec(cvar.u_guide==4);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==5));  dDM.dDM5V(mp.dm5.act_ele) = mp.dm_weights(5)*duVec(cvar.u_guide==5);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==6));  dDM.dDM6V(mp.dm6.act_ele) = mp.dm_weights(6)*duVec(cvar.u_guide==6);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==7));  dDM.dDM7V(mp.dm7.act_ele) = mp.dm_weights(7)*duVec(cvar.u_guide==7);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==8));  dDM.dDM8V(mp.dm8.act_ele) = mp.dm_weights(8)*duVec(cvar.u_guide==8);  end % Parse the command vector to get component for DM and apply the DM's weight
if(any(mp.dm_ind==9));  dDM.dDM9V(mp.dm9.act_ele) = mp.dm_weights(9)*duVec(cvar.u_guide==9);  end % Parse the command vector to get component for DM and apply the DM's weight


 

%% Combine the delta command with the previous command

if(any(mp.dm_ind==1));  mp.dm1.V = cvar.DM1Vnom + dDM.dDM1V;  end
if(any(mp.dm_ind==2));  mp.dm2.V = cvar.DM2Vnom + dDM.dDM2V;  end
if(any(mp.dm_ind==3));  mp.dm3.V = cvar.DM3Vnom + dDM.dDM3V;  end
if(any(mp.dm_ind==4));  mp.dm4.V = cvar.DM4Vnom + dDM.dDM4V;  end
if(any(mp.dm_ind==5));  mp.dm5.V = cvar.DM5Vnom + dDM.dDM5V;  end
if(any(mp.dm_ind==6));  mp.dm6.V = cvar.DM6Vnom + dDM.dDM6V;  end
if(any(mp.dm_ind==7));  mp.dm7.V = cvar.DM7Vnom + dDM.dDM7V;  end
if(any(mp.dm_ind==8));  mp.dm8.V = cvar.DM8Vnom + dDM.dDM8V;  end
if(any(mp.dm_ind==9));  mp.dm9.V = cvar.DM9Vnom + dDM.dDM9V;  end


%% Take images to empirically check contrast at that mu value
Itotal = falco_get_summed_image(mp);
InormAvg = mean(Itotal(mp.F4.corr.maskBool));
cvar.cMin = InormAvg;

cvar.log10regUsed = -10; %--Dummy value



end