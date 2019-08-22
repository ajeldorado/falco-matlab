% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function jac = model_Jacobian_VC(mp, im, whichDM)
%--Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%  This function is for the apodized/unapodized vortex coronagraphs.
%
% REVISION HISTORY:
% --------------
% Modified on 2018-04-22 by A.J. Riggs to use propcustom_relay.
% Modified on 2019-02-12 by G. Ruane to allow scalar vortex models 
% Modified on 2018-03-01 by A.J. Riggs to debug the vortex coronagraph.
% Modified on 2017-11-13 by A.J. Riggs to be compatible with parfor.
% Modified on 2017-11-09 by A.J. Riggs to have the Jacobian calculation be
% its own function.
% Modified on 2017-10-17 by A.J. Riggs to have model_compact.m be a wrapper. All the 
%  actual compact models have been moved to sub-routines for clarity.
% Modified on 19 June 2017 by A.J. Riggs to use lower resolution than the
%   full model.
% Modified by A.J. Riggs on 18 August 2016 from hcil_model.m to model_compact.m.
% Modified by A.J. Riggs on 18 Feb 2015 from HCIL_model_lab_BB_v3.m to hcil_model.m.
% ---------------
%
% INPUTS:
% -mp = structure of model parameters
% -im = index of the pair of sub-bandpass index and Zernike mode index
% -whichDM = which DM number
%
% OUTPUTS:
% -Gzdl = Jacobian for the specified Zernike mode (z), DM (d), and sub-bandpass (l).

function Gzdl = model_Jacobian_VC(mp, im, whichDM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modvar.sbpIndex = mp.jac.sbp_inds(im);
modvar.zernIndex = mp.jac.zern_inds(im);

lambda = mp.sbp_centers(modvar.sbpIndex); 
mirrorFac = 2; % Phase change is twice the DM surface height.f
NdmPad = mp.compact.NdmPad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input E-fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ein = mp.P1.compact.E(:,:,modvar.sbpIndex);  

%--Apply a Zernike (in amplitude) at input pupil
if(modvar.zernIndex~=1)
    indsZnoll = modvar.zernIndex; %--Just send in 1 Zernike mode
    zernMat = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    zernMat = padOrCropEven(zernMat,mp.P1.compact.Narr);
    Ein = Ein.*zernMat*(2*pi/lambda)*mp.jac.Zcoef(mp.jac.zerns==modvar.zernIndex);
    keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pupil = padOrCropEven(mp.P1.compact.mask,NdmPad);
Ein = padOrCropEven(Ein,NdmPad);
%--Re-image the apodizer from pupil P3 back to pupil P2. (Sign of mp.Nrelay2to3 doesn't matter.)
if(mp.flagApod) 
    apodReimaged = padOrCropEven(mp.P3.compact.mask, NdmPad);
    apodReimaged = propcustom_relay(apodReimaged,mp.Nrelay2to3,mp.centering);
else
    apodReimaged = ones(NdmPad); 
end

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.compact.mask, NdmPad); else; DM1stop = ones(NdmPad); end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.compact.mask, NdmPad); else; DM2stop = ones(NdmPad); end

if(any(mp.dm_ind==1)); DM1surf = padOrCropEven(mp.dm1.compact.surfM, NdmPad); else; DM1surf = 0; end 
if(any(mp.dm_ind==2)); DM2surf = padOrCropEven(mp.dm2.compact.surfM, NdmPad); else; DM2surf = 0; end 

if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(mp.dm_ind==1))
        DM1surf = gpuArray(DM1surf);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP2 = propcustom_relay(EP1,mp.Nrelay1to2,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 degrees mp.Nrelay1to2 times.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if(abs(mp.d_P2_dm1)~=0) %--E-field arriving at DM1
    Edm1 = propcustom_PTP(EP2,mp.P2.compact.dx*NdmPad,lambda,mp.d_P2_dm1);
else
    Edm1 = EP2;
end
Edm1 = DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation: 2 DMs, apodizer, binary-amplitude FPM, LS, and final focal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Minimum FPM resolution for Jacobian calculations (in pixels per lambda/D)
minPadFacVortex = 8; 

%--DM1---------------------------------------------------------
if(whichDM==1) 
    if (mp.flagFiber)
        if(mp.flagLenslet)
            Gzdl = zeros(mp.Fend.Nlens,mp.dm1.Nele);
        else
            Gzdl = zeros(mp.Fend.corr.Npix,mp.dm1.Nele);
        end
    else
        Gzdl = zeros(mp.Fend.corr.Npix,mp.dm1.Nele); %--Initialize the Jacobian
    end

    %--Array size for planes P3, F3, and P4
    Nfft1 = 2^(ceil(log2(max([mp.dm1.compact.NdmPad, minPadFacVortex*mp.dm1.compact.Nbox])))); %--Don't crop--but do pad if necessary.
    
    if(numel(mp.F3.VortexCharge)==1)
        charge = mp.F3.VortexCharge;
    else
        charge = interp1(mp.F3.VortexCharge_lambdas,mp.F3.VortexCharge,lambda,'linear','extrap');
    end

    %--Generate vortex FPM with fftshift already applied
    fftshiftVortex = fftshift( falco_gen_vortex_mask( charge, Nfft1) );
    
    %--Two array sizes (at same resolution) of influence functions for MFT and angular spectrum
    NboxPad1AS = mp.dm1.compact.NboxAS; %--Power of 2 array size for FFT-AS propagations from DM1->DM2->DM1
    mp.dm1.compact.xy_box_lowerLeft_AS = mp.dm1.compact.xy_box_lowerLeft - (mp.dm1.compact.NboxAS-mp.dm1.compact.Nbox)/2; %--Adjust the sub-array location of the influence function for the added zero padding

    if(any(mp.dm_ind==2))
        DM2surf = padOrCropEven(DM2surf,mp.dm1.compact.NdmPad);
    else
        DM2surf = zeros(mp.dm1.compact.NdmPad);
    end
    
    if(mp.flagDM2stop)
        DM2stop = padOrCropEven(DM2stop,mp.dm1.compact.NdmPad);
    else
        DM2stop = ones(mp.dm1.compact.NdmPad);
    end
    
    apodReimaged = padOrCropEven( apodReimaged, mp.dm1.compact.NdmPad);
    Edm1pad = padOrCropEven(Edm1,mp.dm1.compact.NdmPad); %--Pad or crop for expected sub-array indexing

    %--Propagate each actuator from DM2 through the optical system
    Gindex = 1; % initialize index counter
    for iact=mp.dm1.act_ele(:).'
        if(any(any(mp.dm1.compact.inf_datacube(:,:,iact))))  %--Only compute for acutators specified for use or for influence functions that are not zeroed out
        
            %--x- and y- coordinates of the padded influence function in the full padded pupil
            x_box_AS_ind = mp.dm1.compact.xy_box_lowerLeft_AS(1,iact):mp.dm1.compact.xy_box_lowerLeft_AS(1,iact)+NboxPad1AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = mp.dm1.compact.xy_box_lowerLeft_AS(2,iact):mp.dm1.compact.xy_box_lowerLeft_AS(2,iact)+NboxPad1AS-1; % y-indices in pupil arrays for the box

            %--Propagate from DM1 to DM2, and then back to P2
            dEbox = (mirrorFac*2*pi*1j/lambda)*padOrCropEven(mp.dm1.VtoH(iact)*mp.dm1.compact.inf_datacube(:,:,iact),NboxPad1AS); %--Pad influence function at DM1 for angular spectrum propagation.
            
            if(mp.useGPU)
                dEbox = gpuArray(dEbox);
            end
            
            dEbox = propcustom_PTP(dEbox.*Edm1pad(y_box_AS_ind,x_box_AS_ind),mp.P2.compact.dx*NboxPad1AS,lambda,mp.d_dm1_dm2); % forward propagate to DM2 and apply DM2 E-field
            dEP2box = propcustom_PTP(dEbox.*DM2stop(y_box_AS_ind,x_box_AS_ind).*exp(mirrorFac*2*pi*1j/lambda*DM2surf(y_box_AS_ind,x_box_AS_ind)),mp.P2.compact.dx*NboxPad1AS,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1) ); % back-propagate to DM1

            %--To simulate going forward to the next pupil plane (with the apodizer) most efficiently, 
            % First, back-propagate the apodizer (by rotating 180-degrees) to the previous pupil.
            dEP2boxEff = apodReimaged(y_box_AS_ind,x_box_AS_ind).*dEP2box; %--Apply 180deg-rotated SP mask.
            
            %--Re-insert the window around the influence function back into the full beam array.
            if(isa(dEP2boxEff,'gpuArray'))
                EP2eff = gpuArray.zeros(mp.dm1.compact.NdmPad);
            else
                EP2eff = zeros(mp.dm1.compact.NdmPad);
            end
            
            EP2eff(y_box_AS_ind,x_box_AS_ind) = dEP2boxEff;
            
            %--Forward propagate from P2 (effective) to P3
            EP3 = propcustom_relay(EP2eff,mp.Nrelay2to3,mp.centering); 
            
            %--Pad pupil P3 for FFT
            EP3pad = padOrCropEven(EP3, Nfft1);

            %--FFT from P3 to Fend.and apply vortex
            EF3 = fftshiftVortex.*fft2(fftshift(EP3pad))/Nfft1;

            %--FFT from Vortex FPM to Lyot Plane
            EP4 = fftshift(fft2(EF3))/Nfft1;
            EP4 = propcustom_relay(EP4,mp.Nrelay3to4-1,mp.centering); %--Add more re-imaging relays if necessary
            if(Nfft1 > mp.P4.compact.Narr)
                EP4 = mp.P4.compact.croppedMask.*padOrCropEven(EP4,mp.P4.compact.Narr); %--Crop EP4 and then apply Lyot stop 
            else
                EP4 = padOrCropEven(mp.P4.compact.croppedMask,Nfft1).*EP4; %--Crop the Lyot stop and then apply it.
            end
            EP4 = propcustom_relay(EP4,mp.NrelayFend,mp.centering); %--Rotate the final image 180 degrees if necessary

            %--MFT to detector
            if(mp.flagFiber)
                if(mp.flagLenslet)
                    for nlens = 1:mp.Fend.Nlens
                        EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering, 'xfc', mp.Fend.x_lenslet_phys(nlens), 'yfc', mp.Fend.y_lenslet_phys(nlens));
                        Elenslet = EFend.*mp.Fend.lenslet.mask;
                        EF5 = propcustom_mft_PtoF(Elenslet, mp.lensletFL, lambda, mp.Fend.dxi, mp.F5.dxi, mp.F5.Nxi, mp.F5.deta, mp.F5.Neta, mp.centering);
                        Gzdl(nlens,Gindex) = max(max(mp.F5.fiberMode(:,:,modvar.sbpIndex))).*sum(sum(mp.F5.fiberMode(:,:,modvar.sbpIndex).*EF5));
                    end
                else
                    EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering);

                    Gzdltemp = zeros(mp.Fend.Nxi, mp.Fend.Neta);
                    for i=1:mp.Fend.Nfiber
                        Eonefiber = mp.Fend.fiberMode(:,:,modvar.sbpIndex,i).*sum(sum(mp.Fend.fiberMode(:,:,modvar.sbpIndex,i).*conj(EFend)));
                        Gzdltemp = Gzdltemp + Eonefiber;
                    end
                    Gzdl(:,Gindex) = Gzdltemp(mp.Fend.corr.inds);
                end
            else    
                EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering);

                if(mp.useGPU)
                    EFend = gather(EFend);
                end
            
                Gzdl(:,Gindex) = EFend(mp.Fend.corr.inds)/sqrt(mp.Fend.compact.I00(modvar.sbpIndex));
            end
        end
        Gindex = Gindex + 1;
    end
end    

%--DM2---------------------------------------------------------
if(whichDM==2)
    if(mp.flagFiber)
        if(mp.flagLenslet)
            Gzdl = zeros(mp.Fend.Nlens,mp.dm2.Nele);
        else
            Gzdl = zeros(mp.Fend.corr.Npix,mp.dm2.Nele);
        end
    else
        Gzdl = zeros(mp.Fend.corr.Npix,mp.dm2.Nele); %--Initialize the Jacobian
    end
    
    %--Array size for planes P3, F3, and P4
    Nfft2 = 2^(ceil(log2(max([mp.dm2.compact.NdmPad, minPadFacVortex*mp.dm2.compact.Nbox])))); %--Don't crop--but do pad if necessary.
    
    if(numel(mp.F3.VortexCharge)==1)
        charge = mp.F3.VortexCharge;
    else
        charge = interp1(mp.F3.VortexCharge_lambdas,mp.F3.VortexCharge,lambda,'linear','extrap');
    end

    %--Generate vortex FPM with fftshift already applied
    fftshiftVortex = fftshift( falco_gen_vortex_mask( charge, Nfft2) );
    
    %--Two array sizes (at same resolution) of influence functions for MFT and angular spectrum
    NboxPad2AS = mp.dm2.compact.NboxAS; 
    mp.dm2.compact.xy_box_lowerLeft_AS = mp.dm2.compact.xy_box_lowerLeft - (NboxPad2AS-mp.dm2.compact.Nbox)/2; %--Account for the padding of the influence function boxes
    
    apodReimaged = padOrCropEven( apodReimaged, mp.dm2.compact.NdmPad);
    DM2stop = padOrCropEven(DM2stop,mp.dm2.compact.NdmPad);
        
    %--Propagate full field to DM2 before back-propagating in small boxes
    Edm2inc = padOrCropEven( propcustom_PTP(Edm1,mp.compact.NdmPad*mp.P2.compact.dx,lambda,mp.d_dm1_dm2), mp.dm2.compact.NdmPad); % E-field incident upon DM2
    Edm2inc = padOrCropEven(Edm2inc,mp.dm2.compact.NdmPad);
    Edm2 = DM2stop.*Edm2inc.*exp(mirrorFac*2*pi*1j/lambda*padOrCropEven(DM2surf,mp.dm2.compact.NdmPad)); % Initial E-field at DM2 including its own phase contribution
    
    %--Propagate each actuator from DM2 through the rest of the optical system
    Gindex = 1; % initialize index counter
    for iact=mp.dm2.act_ele(:).'
        if(any(any(mp.dm2.compact.inf_datacube(:,:,iact))) ) 
            
            %--x- and y- coordinates of the padded influence function in the full padded pupil
            x_box_AS_ind = mp.dm2.compact.xy_box_lowerLeft_AS(1,iact):mp.dm2.compact.xy_box_lowerLeft_AS(1,iact)+NboxPad2AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = mp.dm2.compact.xy_box_lowerLeft_AS(2,iact):mp.dm2.compact.xy_box_lowerLeft_AS(2,iact)+NboxPad2AS-1; % y-indices in pupil arrays for the box

            dEbox = mp.dm2.VtoH(iact)*(mirrorFac*2*pi*1j/lambda)*padOrCropEven(mp.dm2.compact.inf_datacube(:,:,iact),NboxPad2AS); %--the padded influence function at DM2
            
            if(mp.useGPU)
                dEbox = gpuArray(dEbox);
            end
            
            dEP2box = propcustom_PTP(dEbox.*Edm2(y_box_AS_ind,x_box_AS_ind),mp.P2.compact.dx*NboxPad2AS,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); % back-propagate to pupil P2

            %--To simulate going forward to the next pupil plane (with the apodizer) most efficiently, 
            % First, back-propagate the apodizer (by rotating 180-degrees) to the previous pupil.
            dEP2boxEff = apodReimaged(y_box_AS_ind,x_box_AS_ind).*dEP2box; %--Apply de-rotated SP mask.

            if(isa(dEP2boxEff,'gpuArray'))
                EP2eff = gpuArray.zeros(mp.dm2.compact.NdmPad);
            else
                EP2eff = zeros(mp.dm2.compact.NdmPad);
            end
            
            EP2eff(y_box_AS_ind,x_box_AS_ind) = dEP2boxEff;

            %--Forward propagate from P2 (effective) to P3
            EP3 = propcustom_relay(EP2eff,mp.Nrelay2to3,mp.centering); 

            %--Pad pupil P3 for FFT
            EP3pad = padOrCropEven(EP3, Nfft2);

            %--FFT from P3 to Fend.and apply vortex
            EF3 = fftshiftVortex.*fft2(fftshift(EP3pad))/Nfft2;

            %--FFT from Vortex FPM to Lyot Plane
            EP4 = fftshift(fft2(EF3))/Nfft2;
            EP4 = propcustom_relay(EP4,mp.Nrelay3to4-1,mp.centering); %--Add more re-imaging relays if necessary
            if(Nfft2 > mp.P4.compact.Narr)
                EP4 = mp.P4.compact.croppedMask.*padOrCropEven(EP4,mp.P4.compact.Narr); %--Crop EP4 and then apply Lyot stop 
            else
                EP4 = padOrCropEven(mp.P4.compact.croppedMask,Nfft2).*EP4; %--Crop the Lyot stop and then apply it.
            end
            EP4 = propcustom_relay(EP4,mp.NrelayFend,mp.centering); %--Rotate the final image 180 degrees if necessary

            %--MFT to detector
            if(mp.flagFiber)
                if(mp.flagLenslet)
                    for nlens = 1:mp.Fend.Nlens
                        EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering, 'xfc', mp.Fend.x_lenslet_phys(nlens), 'yfc', mp.Fend.y_lenslet_phys(nlens));
                        Elenslet = EFend.*mp.Fend.lenslet.mask;
                        EF5 = propcustom_mft_PtoF(Elenslet, mp.lensletFL, lambda, mp.Fend.dxi, mp.F5.dxi, mp.F5.Nxi, mp.F5.deta, mp.F5.Neta, mp.centering);
                        Gzdl(nlens,Gindex) = max(max(mp.F5.fiberMode(:,:,modvar.sbpIndex))).*sum(sum(mp.F5.fiberMode(:,:,modvar.sbpIndex).*EF5));
                    end
                else
                    EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering);
                    Gzdltemp = zeros(mp.Fend.Nxi, mp.Fend.Neta);
                    for i=1:mp.Fend.Nfiber
                        Eonefiber = mp.Fend.fiberMode(:,:,modvar.sbpIndex,i).*sum(sum(mp.Fend.fiberMode(:,:,modvar.sbpIndex,i).*conj(EFend)));
                        Gzdltemp = Gzdltemp + Eonefiber;
                    end
                    Gzdl(:,Gindex) = Gzdltemp(mp.Fend.corr.inds);
                end
            else    
                EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering);

                if(mp.useGPU)
                    EFend = gather(EFend);
                end
            
                Gzdl(:,Gindex) = EFend(mp.Fend.corr.inds)/sqrt(mp.Fend.compact.I00(modvar.sbpIndex));
            end
        end
        Gindex = Gindex + 1;
    end
end

end % End of function