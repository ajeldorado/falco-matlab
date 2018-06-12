% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_compact_LC(mp, DM, modvar)
%--Blind model used by the estimator and controller
%  Does not include unknown aberrations/errors that are in the full model.
%
%  This function is for the Lyot coronagraph, DMLC, and APLC.
%
% REVISION HISTORY:
% --------------
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

% function Eout = model_compact_LC(mp, DM, modvar)
% 

function Eout = model_compact_LC(mp, DM, modvar)
lambda = mp.sbp_center_vec(modvar.sbpIndex);
mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = DM.compact.NdmPad;

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
if(any(DM.dm_ind==1)); DM1surf = falco_gen_dm_surf(DM.dm1, DM.dm1.compact.dx, NdmPad); else DM1surf = 0; end %--Pre-compute the starting DM1 surface
if(any(DM.dm_ind==2)); DM2surf = falco_gen_dm_surf(DM.dm2, DM.dm2.compact.dx, NdmPad); else DM2surf = 0; end %--Pre-compute the starting DM2 surface

pupil = padOrCropEven(mp.P1.compact.mask,NdmPad);
Ein = padOrCropEven(Ein,DM.compact.NdmPad);

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.compact.mask, NdmPad); else DM1stop = 1; end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.compact.mask, NdmPad); else DM2stop = 1; end

if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(DM.dm_ind==1)); DM1surf = gpuArray(DM1surf); end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation: 2 DMs, apodizer, binary-amplitude FPM, LS, and final focal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP2 = propcustom_2FT(EP1,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 deg.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if( abs(mp.d_P2_dm1)~=0 ); Edm1 = propcustom_PTP(EP2,mp.P2.compact.dx*NdmPad,lambda,mp.d_P2_dm1); else Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.compact.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if( mp.d_P2_dm1 + mp.d_dm1_dm2 == 0 ); EP2eff = Edm2; else EP2eff = propcustom_PTP(Edm2,mp.P2.compact.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); end %--Back propagate to pupil P2

%--Rotate 180 degrees to propagate to pupil P3
EP3 = propcustom_2FT(EP2eff, mp.centering);

%--Apply apodizer mask.
if(mp.flagApod)
    EP3 = mp.P3.compact.mask.*padOrCropEven(EP3, mp.P3.compact.Narr); 
end

%--MFT from SP to FPM (i.e., P3 to F3)
EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.compact.dx,mp.F3.compact.dxi,mp.F3.compact.Nxi,mp.F3.compact.deta,mp.F3.compact.Neta,mp.centering); %--E-field incident upon the FPM
EF3 = (1-mp.F3.compact.mask.amp).*EF3inc; %--Apply (1-FPM) for Babinet's principle later

%--MFT from FPM to Lyot Plane (i.e., F3 to P4)
EP4sub = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.compact.dxi,mp.F3.compact.deta,mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering); % Subtrahend term for Babinet's principle     

%--Use Babinet's principle at the Lyot plane.
EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.compact.Narr); %--Crop down to the size of the Lyot stop opening
EP4 = mp.P4.compact.croppedMask.*(EP4noFPM-EP4sub); % (1j)^2 comes from the coefficients of the 2 skipped MFTs to reimage this pupil

if(isfield(modvar,'flagGetNormVal'))
    if(modvar.flagGetNormVal==true)
        EP4 = mp.P4.compact.croppedMask.*(EP4noFPM);  %--No FPM, so include just the effect of the DMs and LS
    end
end

% DFT to camera
EF4 = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.F4.compact.dxi,mp.F4.compact.Nxi,mp.F4.compact.deta,mp.F4.compact.Neta,mp.centering);


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


    
