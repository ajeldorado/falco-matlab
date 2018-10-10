% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_full_SPHLC(mp,   modvar)
%--Full-knowledge optical model.
%    --> Not used by the estimator and controller.
%    --> Only used to create simulated intensity images.
%
% New or changed variables in going from HLC to SPHLC:
%  mp.F3.full.dxi  --> mp.F3.full.in.dxi,  mp.F3.full.out.dxi
%  mp.F3.full.deta --> mp.F3.full.in.deta, mp.F3.full.out.deta
%  mp.F3.full.Nxi  --> mp.F3.full.in.Nxi,  mp.F3.full.out.Nxi
%  mp.F3.full.Neta --> mp.F3.full.in.Neta, mp.F3.full.out.Neta
%  FPM --> FPMinner
%  mp.F3.full.iris (new)
%
% 
%
% REVISION HISTORY:
% --------------
% Modified on 2018-05-29 by A.J. Riggs from the HLC model to the SPHLC
% model. 
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


function Eout = model_full_SPHLC(mp,   modvar)

lambda = mp.sbp_centers(modvar.sbpIndex)*mp.lamFac_vec(modvar.wpsbpIndex);
mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.full.NdmPad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(any(mp.dm_ind==1)) 
    if( isfield(mp.dm1,'surfM') )
        DM1surf = padOrCropEven(mp.dm1.surfM, NdmPad);
    else
        DM1surf = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,NdmPad); 
    end
else
    DM1surf = 0;
end

if(any(mp.dm_ind==2))
    if( isfield(mp.dm2,'surfM') )   
        DM2surf = padOrCropEven(mp.dm2.surfM, NdmPad);
    else
        DM2surf = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,NdmPad); 
    end
else
    DM2surf = 0;
end

%--Complex transmission map of the FPM.
ilam = (modvar.sbpIndex-1)*mp.Nwpsbp + modvar.wpsbpIndex;
if( isfield(mp,'FPMcubeFull') )  %--Load it if stored
    FPMinner = mp.FPMcubeFull(:,:,ilam); %padOrCropEven(mp.FPMcubeFull, mp.dm9.NxiFPM);
else %--Otherwise generate it
    FPMinner = falco_gen_HLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'full'); %padOrCropEven( ,mp.dm9.NxiFPM);
end

%--Complex transmission of the points outside the FPM (just fused silica with neither dielectric nor metal).
ind_metal = falco_discretize_FPM_surf(0, mp.t_metal_nm_vec, mp.dt_metal_nm); %--Obtain the indices of the nearest thickness values in the complex transmission datacube.
ind_diel = falco_discretize_FPM_surf(0, mp.t_diel_nm_vec,  mp.dt_diel_nm); %--Obtain the indices of the nearest thickness values in the complex transmission datacube.
transOuterFPM = mp.complexTransFull(ind_diel,ind_metal,ilam); %--Complex transmission of the points outside the FPM (just fused silica with neither dielectric nor metal).            

% %--DEBUGGING: VERIFYING THAT THE FPM PHASE DOES NOT SHIFT ARBITRARILY
% figure(411); imagesc(abs(FPM)); axis xy equal tight; colorbar;
% figure(412); imagesc(angle(FPM)); axis xy equal tight; colorbar;
% DMcopy = DM;
% DMcopy.dm9.V(round(mp.dm9.Nele/2)+30) = 400;
% FPM2 = falco_gen_HLC_FPM_complex_trans_mat(DMcopy,mp,modvar.sbpIndex,modvar.wpsbpIndex,'full'); %padOrCropEven( ,mp.dm9.NxiFPM);
% figure(511); imagesc(abs(FPM2)); axis xy equal tight; colorbar;
% figure(512); imagesc(angle(FPM2)); axis xy equal tight; colorbar;
% figure(513); imagesc(abs(FPM-FPM2)); axis xy equal tight; colorbar;
% figure(514); imagesc(angle(FPM-FPM2)); axis xy equal tight; colorbar;


pupil = padOrCropEven(mp.P1.full.mask,NdmPad);
Ein = padOrCropEven(Ein,NdmPad);

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.full.mask, NdmPad); else; DM1stop = 1; end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.full.mask, NdmPad); else; DM2stop = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation: entrance pupil, 2 DMs, (optional) apodizer, complex-valued FPM, LS, and final focal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP1 = EP1.*padOrCropEven(mp.P3.full.mask,length(EP1)); mp.flagApod = false;%--TESTING: put apodizer before the DMs
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



if(mp.flagSPHLConeFPM) %--Only inner region
    
    %--MFT from SP to FPM (i.e., P3 to F3)
    EF3inInc  = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx, mp.F3.full.in.dxi, mp.F3.full.in.Nxi, mp.F3.full.in.deta, mp.F3.full.in.Neta,  mp.centering);

    % Apply FPM (inner region is complex; iris for the outer part)
    EF3in = mp.F3.full.irisHD.*FPMinner.*EF3inInc; %-transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.

    %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
    EP4 = propcustom_mft_FtoP(EF3in, mp.fl,lambda, mp.F3.full.in.dxi,  mp.F3.full.in.deta,  mp.P4.full.dx,mp.P4.full.Narr,mp.centering);

    %--Do NOT apply (inner) FPM if normalization value is being found
    if(isfield(modvar,'flagGetNormVal'))
        if(modvar.flagGetNormVal==true)
            %--MFT from SP to FPM (i.e., P3 to F3)
            EF3outInc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx, mp.F3.full.out.dxi,mp.F3.full.out.Nxi,mp.F3.full.out.deta,mp.F3.full.out.Neta, mp.centering);
            % Apply FPM (inner region is complex; iris for the outer part)
            EF3out = transOuterFPM*mp.F3.full.iris.*EF3outInc;
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4 = propcustom_mft_FtoP(EF3out,mp.fl,lambda, mp.F3.full.out.dxi, mp.F3.full.out.deta, mp.P4.full.dx,mp.P4.full.Narr,mp.centering);
        end
    end
    
else

    %--MFT from apodizer plane to FPM (i.e., P3 to F3) in two parts
    EF3inInc  = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx, mp.F3.full.in.dxi, mp.F3.full.in.Nxi, mp.F3.full.in.deta, mp.F3.full.in.Neta,  mp.centering);
    EF3outInc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx, mp.F3.full.out.dxi,mp.F3.full.out.Nxi,mp.F3.full.out.deta,mp.F3.full.out.Neta, mp.centering);

    % Apply FPM (inner region is complex; iris for the outer part)
    EF3in = (FPMinner - transOuterFPM).*EF3inInc; %-transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
    EF3out = transOuterFPM*mp.F3.full.iris.*EF3outInc;

    %--Do NOT apply (inner) FPM if normalization value is being found
    if(isfield(modvar,'flagGetNormVal'))
        if(modvar.flagGetNormVal==true)
            EF3in = 0;
        end
    end

    %--MFT from FPM to Lyot Plane (i.e., F3 to P4) in 2 parts
    EP4 = propcustom_mft_FtoP(EF3in, mp.fl,lambda, mp.F3.full.in.dxi,  mp.F3.full.in.deta,  mp.P4.full.dx,mp.P4.full.Narr,mp.centering)+...
          propcustom_mft_FtoP(EF3out,mp.fl,lambda, mp.F3.full.out.dxi, mp.F3.full.out.deta, mp.P4.full.dx,mp.P4.full.Narr,mp.centering);

end


% %--Apply Lyot stop
% IP4 = abs(EP4).^2;
% IP4 = IP4/max(max(IP4));
% figure; colorbar; imagesc(log10(IP4),[-4 0]); colorbar;
% EP4 = mp.P4.full.croppedMask.*EP4; 
% figure; colorbar; imagesc(log10(abs(EP4).^2/max(max(IP4))),[-4 0]); colorbar;
% 
% keyboard


%--Apply Lyot stop
EP4 = mp.P4.full.croppedMask.*EP4; 


%--MFT from Lyot Stop to final focal plane (i.e., P4 to F4)
EF4 = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.full.dx,mp.F4.dxi,mp.F4.Nxi,mp.F4.deta,mp.F4.Neta,mp.centering);

%--Don't apply FPM if normalization value is being found, or if the flag doesn't exist (for testing only)
Eout = EF4; %--Don't normalize if normalization value is being found
if(isfield(modvar,'flagGetNormVal'))
    if(modvar.flagGetNormVal==false)
        Eout = EF4/sqrt(mp.F4.full.I00(modvar.sbpIndex)); %--Apply normalization
    end
elseif(isfield(mp.F4.full,'I00'))
    Eout = EF4/sqrt(mp.F4.full.I00(modvar.sbpIndex)); %--Apply normalization
end


end %--END OF FUNCTION


    
