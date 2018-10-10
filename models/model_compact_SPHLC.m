% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_compact_SPHLC(mp, DM, modvar)
%--Blind model used by the estimator and controller
%  Does not include unknown aberrations/errors that are in the full model.
%
%
% New or changed variables in going from HLC to SPHLC:
%  mp.F3.compact.dxi  --> mp.F3.compact.in.dxi,  mp.F3.compact.out.dxi
%  mp.F3.compact.deta --> mp.F3.compact.in.deta, mp.F3.compact.out.deta
%  mp.F3.compact.Nxi  --> mp.F3.compact.in.Nxi,  mp.F3.compact.out.Nxi
%  mp.F3.compact.Neta --> mp.F3.compact.in.Neta, mp.F3.compact.out.Neta
%  FPM --> FPMinner
%  mp.F3.compact.iris (new)
%
%
% REVISION HISTORY:
% --------------
% Modified on 2018-05-29 by A.J. Riggs from the HLC model to the SPHLC
% model.
% Modified on 2018-01-23 by A.J. Riggs to allow DM1 to not be at a pupil
%  and to have an aperture stop.
% Modified on 2017-11-09 by A.J. Riggs to remove the Jacobian calculation.
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
% -Eout = electric field in the final focal plane
%
% modvar structure fields (4):
% -sbpIndex
% -wpsbpIndex
% -whichSource
% -flagGenMat

% function Eout = model_compact_HLC(mp, DM, modvar)
% 

function Eout = model_compact_SPHLC(mp, DM, modvar)

lambda = mp.sbp_centers(modvar.sbpIndex); 
mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.compact.NdmPad;

if(flagEval)
    dxi = mp.F4.eval.dxi;
    Nxi = mp.F4.eval.Nxi;
    deta = mp.F4.eval.deta;
    Neta = mp.F4.eval.Neta; 
else
    dxi = mp.F4.dxi;
    Nxi = mp.F4.Nxi;
    deta = mp.F4.deta;
    Neta = mp.F4.Neta; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input E-fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Include the tip/tilt in the input wavefront
if(isfield(mp,'ttx'))
    %--Scale by lambda/lambda0 because ttx and tty are in lambda0/D
    x_offset = mp.ttx(modvar.ttIndex);
    y_offset = mp.tty(modvar.ttIndex);

    TTphase = (-1)*(2*pi*(x_offset*mp.P2.compact.XsDL + y_offset*mp.P2.compact.YsDL));
    Ett = exp(1i*TTphase*mp.lambda0/lambda);
    Ein = Ett.*mp.P1.compact.E(:,:,modvar.sbpIndex);  

else %--Backward compatible with code without tip/tilt offsets in the Jacobian
    Ein = mp.P1.compact.E(:,:,modvar.sbpIndex);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Compute the DM surfaces for the current DM commands
if(any(mp.dm_ind==1)); DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, NdmPad); else; DM1surf = 0; end %--Pre-compute the starting DM1 surface
if(any(mp.dm_ind==2)); DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, NdmPad); else; DM2surf = 0; end %--Pre-compute the starting DM2 surface
% if(any(mp.dm_ind==9)); DM9phase = padOrCropEven(falco_dm_surf_from_cube(mp.dm9,mp.dm9.compact),mp.F3.compact.Nxi); else DM9phase = 0; end %--Pre-compute the starting DM9 surface
   
%--Generate the FPM's complex transmission map
FPMinner = falco_gen_HLC_FPM_complex_trans_mat(DM,mp,modvar.sbpIndex,modvar.wpsbpIndex,'compact'); 
% figure(421); imagesc(abs(FPM)); axis xy equal tight; colorbar;
% figure(422); imagesc(angle(FPM)); axis xy equal tight; colorbar;
% FPM = padOrCropEven( FPM,mp.dm9.compact.NxiFPM);

%--Complex transmission of the points outside the FPM (just fused silica with neither dielectric nor metal).
ilam = modvar.sbpIndex;
ind_metal = falco_discretize_FPM_surf(0, mp.t_metal_nm_vec, mp.dt_metal_nm); %--Obtain the indices of the nearest thickness values in the complex transmission datacube.
ind_diel = falco_discretize_FPM_surf(0, mp.t_diel_nm_vec,  mp.dt_diel_nm); %--Obtain the indices of the nearest thickness values in the complex transmission datacube.
transOuterFPM = mp.complexTransCompact(ind_diel,ind_metal,ilam); %--Complex transmission of the points outside the FPM (just fused silica with neither dielectric nor metal).            

% %--DEBUGGING: VERIFYING THAT THE FPM PHASE DOES NOT SHIFT ARBITRARILY
% figure(411); imagesc(abs(FPM)); axis xy equal tight; colorbar;
% figure(412); imagesc(angle(FPM)); axis xy equal tight; colorbar;
% DMcopy = DM;
% DMcopy.dm9.V(round(mp.dm9.Nele/2)+30) = 400;
% FPM2 = falco_gen_HLC_FPM_complex_trans_mat(DMcopy,mp,modvar.sbpIndex,modvar.wpsbpIndex,'compact'); 
% figure(511); imagesc(abs(FPM2)); axis xy equal tight; colorbar;
% figure(512); imagesc(angle(FPM2)); axis xy equal tight; colorbar;
% figure(513); imagesc(abs(FPM-FPM2)); axis xy equal tight; colorbar;
% figure(514); imagesc(angle(FPM-FPM2)); axis xy equal tight; colorbar;





pupil = padOrCropEven(mp.P1.compact.mask,NdmPad);
Ein = padOrCropEven(Ein,mp.compact.NdmPad);

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.compact.mask, NdmPad); else; DM1stop = 1; end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.compact.mask, NdmPad); else; DM2stop = 1; end

if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(mp.dm_ind==1)); DM1surf = gpuArray(DM1surf); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation: 2 DMs, apodizer, filtered complex FPM, LS, and final focal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP1 = EP1.*padOrCropEven(mp.P3.compact.mask,length(EP1)); mp.flagApod = false; %--TESTING: put apodizer before the DMs
EP2 = propcustom_2FT(EP1,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 deg.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if( abs(mp.d_P2_dm1)~=0 ); Edm1 = propcustom_PTP(EP2,mp.P2.compact.dx*NdmPad,lambda,mp.d_P2_dm1); else; Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.compact.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if( mp.d_P2_dm1 + mp.d_dm1_dm2 == 0 ); EP2eff = Edm2; else; EP2eff = propcustom_PTP(Edm2,mp.P2.compact.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); end %--Back propagate to pupil P2

%--Rotate 180 degrees to propagate to pupil P3
EP3 = propcustom_2FT(EP2eff, mp.centering);

%--Apply apodizer or pupil mask (if there is one)
if(mp.flagApod)
    EP3 = mp.P3.compact.mask.*padOrCropEven(EP3, mp.P3.compact.Narr); 
end



if(mp.flagSPHLConeFPM) %--Only inner region
    
    %--MFT from SP to FPM (i.e., P3 to F3)
    EF3inInc  = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.compact.dx, mp.F3.compact.in.dxi, mp.F3.compact.in.Nxi, mp.F3.compact.in.deta, mp.F3.compact.in.Neta,  mp.centering);

    % Apply FPM (inner region is complex; iris for the outer part)
    EF3in = mp.F3.compact.irisHD.*FPMinner.*EF3inInc; %-transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.

    %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
    EP4 = propcustom_mft_FtoP(EF3in, mp.fl,lambda, mp.F3.compact.in.dxi,  mp.F3.compact.in.deta,  mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering);

    %--Do NOT apply (inner) FPM if normalization value is being found
    if(isfield(modvar,'flagGetNormVal'))
        if(modvar.flagGetNormVal==true)
            %--MFT from SP to FPM (i.e., P3 to F3)
            EF3outInc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.compact.dx, mp.F3.compact.out.dxi,mp.F3.compact.out.Nxi,mp.F3.compact.out.deta,mp.F3.compact.out.Neta, mp.centering);
            % Apply FPM (inner region is complex; iris for the outer part)
            EF3out = transOuterFPM*mp.F3.compact.iris.*EF3outInc;
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4 = propcustom_mft_FtoP(EF3out,mp.fl,lambda, mp.F3.compact.out.dxi, mp.F3.compact.out.deta, mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering);
        end
    end

    
else %--FPM computed as two separate parts: inner region is complex, and iris for the outer part.
    
    %--MFT from SP to FPM (i.e., P3 to F3)
    EF3inInc  = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.compact.dx, mp.F3.compact.in.dxi, mp.F3.compact.in.Nxi, mp.F3.compact.in.deta, mp.F3.compact.in.Neta,  mp.centering);
    EF3outInc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.compact.dx, mp.F3.compact.out.dxi,mp.F3.compact.out.Nxi,mp.F3.compact.out.deta,mp.F3.compact.out.Neta, mp.centering);

    % Apply FPM (inner region is complex; iris for the outer part)
    EF3in = (FPMinner - transOuterFPM).*EF3inInc; %-transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
    EF3out = transOuterFPM*mp.F3.compact.iris.*EF3outInc;

    %--Do NOT apply (inner) FPM if normalization value is being found
    if(isfield(modvar,'flagGetNormVal'))
        if(modvar.flagGetNormVal==true)
            EF3in = 0;
        end
    end

    %--MFT from FPM to Lyot Plane (i.e., F3 to P4) in 2 parts
    EP4 = propcustom_mft_FtoP(EF3in, mp.fl,lambda, mp.F3.compact.in.dxi,  mp.F3.compact.in.deta,  mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering)+...
          propcustom_mft_FtoP(EF3out,mp.fl,lambda, mp.F3.compact.out.dxi, mp.F3.compact.out.deta, mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering);

end


%--Apply Lyot Stop
EP4 = mp.P4.compact.croppedMask.*EP4; %--Apply Lyot stop      

% DFT to camera
EF4 = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.F4.dxi,mp.F4.Nxi,mp.F4.deta,mp.F4.Neta,mp.centering);


%--Don't apply FPM if normalization value is being found, or if the flag doesn't exist (for testing only)
Eout = EF4; %--Don't normalize if normalization value is being found
if(isfield(modvar,'flagGetNormVal'))
    if(modvar.flagGetNormVal==false)
        Eout = EF4/sqrt(mp.F4.compact.I00(modvar.sbpIndex)); %--Apply normalization
    end
elseif(isfield(mp.F4.compact,'I00'))
    Eout = EF4/sqrt(mp.F4.compact.I00(modvar.sbpIndex)); %--Apply normalization
end

if(mp.useGPU)
    Eout = gather(Eout);
end


end %--END OF FUNCTION


    
