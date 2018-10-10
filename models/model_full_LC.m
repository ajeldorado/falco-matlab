% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_full_LC(mp,   modvar)
%--Full-knowledge optical model.
%    --> Not used by the estimator and controller.
%    --> Only used to create simulated intensity images.
%
%  This function is for the Lyot coronagraph, DMLC, and APLC.
%
% REVISION HISTORY:
% --------------
% Modified on 2018-01-23 by A.J. Riggs to allow DM1 to not be at a pupil
%  and to have an aperture stop.
% Modified on 2017-10-17 by A.J. Riggs to have model_full.m be a wrapper. All the 
%  actual full models, including this one, have been moved to sub-routines for clarity.
% Modified by A.J. Riggs from hcil_simTestbed.m to model_full.m.
% Modified on 2015-02-18 by A.J. Riggs from hcil_model.m to hcil_simTestbed.m to inclue
%  extra errors in the model to simulate the actual testbed for generated images.
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


function Eout = model_full_LC(mp,   lambda, Ein, normFac)
% function Eout = model_full_LC(mp,   modvar)

% lambda = mp.sbp_centers(modvar.sbpIndex)*mp.lamFac_vec(modvar.wpsbpIndex);
mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.full.NdmPad;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(any(mp.dm_ind==1)) 
    if( isfield(mp.dm1,'surfM') );   DM1surf = padOrCropEven(mp.dm1.surfM, NdmPad);
    else                            DM1surf = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,NdmPad); end
else DM1surf = 0;
end
if(any(mp.dm_ind==2)) 
    if( isfield(mp.dm2,'surfM') );   DM2surf = padOrCropEven(mp.dm2.surfM, NdmPad);
    else                            DM2surf = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,NdmPad); end
else DM2surf = 0;
end

pupil = padOrCropEven(mp.P1.full.mask,NdmPad);
Ein = padOrCropEven(Ein,NdmPad);

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.full.mask, NdmPad); else; DM1stop = 1; end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.full.mask, NdmPad); else; DM2stop = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation: entrance pupil, 2 DMs, (optional) apodizer, complex-valued FPM, LS, and final focal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mp.dm1.V = 100*(eye(mp.dm1.Nact)+rot90(eye(mp.dm1.Nact),1)); %--DEBUGGING ONLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% DM1surf = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,NdmPad);%--DEBUGGING ONLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP2 = propcustom_2FT(EP1,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 deg.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if( abs(mp.d_P2_dm1)~=0 ); Edm1 = propcustom_PTP(EP2,mp.P2.full.dx*NdmPad,lambda,mp.d_P2_dm1); else; Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.full.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if( mp.d_P2_dm1 + mp.d_dm1_dm2 == 0 ); EP2eff = Edm2; else; EP2eff = propcustom_PTP(Edm2,mp.P2.full.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); end %--Back propagate to pupil P2

%--Rotate 180 degrees to propagate to pupil P3
EP3 = propcustom_2FT(EP2eff, mp.centering);

%--Apply the apodizer mask
if(mp.flagApod)
    EP3 = padOrCropEven(mp.P3.full.mask,mp.P3.full.Narr).*padOrCropEven(EP3, mp.P3.full.Narr); %--Apply apodizer mask.
end

%--MFT from apodizer plane to FPM (i.e., P3 to F3)
EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);

% Apply (1-FPM) for Babinet's principle later
EF3 = (1-mp.F3.full.mask.amp).*EF3inc;

% Use Babinet's principle at the Lyot plane. This is the term without the FPM.
EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.full.Narr);

if(normFac==0) %--Do NOT apply FPM if normalization value is being found
    EP4 = mp.P4.full.croppedMask.*EP4noFPM;
else %--Otherwise apply FPM
    %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
    EP4subtra = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
    %--Babinet's principle at P4
    EP4 = mp.P4.full.croppedMask.*(EP4noFPM-EP4subtra);
end

% %--Do NOT apply FPM if normalization value is being found
% if(isfield(modvar,'flagGetNormVal'))
%     if(modvar.flagGetNormVal==true)
%         EP4 = mp.P4.full.croppedMask.*EP4noFPM;
%     end
% end    

%--MFT from Lyot Stop to final focal plane (i.e., P4 to F4)
EF4 = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.full.dx,mp.F4.dxi,mp.F4.Nxi,mp.F4.deta,mp.F4.Neta,mp.centering);

%--Don't apply FPM if normalization value is being found
if(normFac==0)
    Eout = EF4; %--Don't normalize if normalization value is being found
else
    Eout = EF4/sqrt(normFac); %--Apply normalization
end

% %--Don't apply FPM if normalization value is being found, or if the flag doesn't exist (for testing only)
% Eout = EF4; %--Don't normalize if normalization value is being found
% if(isfield(modvar,'flagGetNormVal'))
%     if(modvar.flagGetNormVal==false)
%         Eout = EF4/sqrt(mp.F4.full.I00(modvar.sbpIndex)); %--Apply normalization
%     end
% elseif(isfield(mp.F4.full,'I00'))
%     Eout = EF4/sqrt(mp.F4.full.I00(modvar.sbpIndex)); %--Apply normalization
% end


end %--END OF FUNCTION


    
