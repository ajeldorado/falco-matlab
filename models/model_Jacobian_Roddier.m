% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function jac = model_Jacobian_Roddier(mp, DM, tsi, whichDM)
%--Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%  This function is for a Roddier coronagraph.
%
% REVISION HISTORY:
% --------------
% Modified on 2017-11-20 from model_Jacobian_LC by G. Ruane
% Modified on 2017-11-13 by A.J. Riggs to be compatible with parfor.
% Modified on 2017-11-09 by A.J. Riggs to compute only one row of the whole Jacobian. 
%  This enables much easier parallelization.
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
% -DM = structure of DM settings
% -tsi = index of the pair of sub-bandpass index and tip/tilt offset index
% -whichDM = which DM number
%
% OUTPUTS:
% -Gttlam = Jacobian for the specified DM and specified T/T-wavelength pair
%

function Gzdl = model_Jacobian_Roddier(mp, im, whichDM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modvar.sbpIndex = mp.jac.sbp_inds(im);
modvar.zernIndex = mp.jac.zern_inds(im);
indZernVec = find(mp.jac.zerns==modvar.zernIndex);

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
    % figure(1); imagesc(zernMat); axis xy equal tight; colorbar; 
    Ein = Ein.*zernMat*(2*pi*1i/lambda)*mp.jac.Zcoef(indZernVec);
end
% figure(71); imagesc(real(Ein)); axis xy equal tight; colorbar; 
% figure(72); imagesc(imag(Ein)); axis xy equal tight; colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pupil = padOrCropEven(mp.P1.compact.mask,NdmPad);
Ein = padOrCropEven(Ein,NdmPad);
if(mp.flagApod) 
    apodRot180 = padOrCropEven( rot90(mp.P3.compact.mask,2), NdmPad );
    if( strcmpi(mp.centering,'pixel') ); apodRot180 = circshift(apodRot180,[1 1]); end %--To undo center offset when pixel centered and rotating by 180 degrees.
else
    apodRot180 = ones(NdmPad); 
end


if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.compact.mask, NdmPad); else; DM1stop = ones(NdmPad); end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.compact.mask, NdmPad); else; DM2stop = ones(NdmPad); end

if(any(mp.dm_ind==1)); DM1surf = padOrCropEven(mp.dm1.compact.surfM, NdmPad);  else; DM1surf = 0; end 
if(any(mp.dm_ind==2)); DM2surf = padOrCropEven(mp.dm2.compact.surfM, NdmPad);  else; DM2surf = 0; end 

if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(mp.dm_ind==1)); DM1surf = gpuArray(DM1surf); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP2 = propcustom_2FT(EP1,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 deg.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if( abs(mp.d_P2_dm1)~=0 ); Edm1 = propcustom_PTP(EP2,mp.P2.compact.dx*NdmPad,lambda,mp.d_P2_dm1); else; Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--DM1---------------------------------------------------------
if(whichDM==1) 
    Gzdl = zeros(mp.F4.corr.Npix,mp.dm1.Nele);
    
    %--Two array sizes (at same resolution) of influence functions for MFT and angular spectrum
    Nbox1 = mp.dm1.compact.Nbox; %--Smaller array size for MFT to FPM after FFT-AS propagations from DM1->DM2->DM1
    NboxPad1AS = mp.dm1.compact.NboxAS; %NboxPad1;%2.^ceil(log2(NboxPad1)); %--Power of 2 array size for FFT-AS propagations from DM1->DM2->DM1
    mp.dm1.compact.xy_box_lowerLeft_AS = mp.dm1.compact.xy_box_lowerLeft - (mp.dm1.compact.NboxAS-mp.dm1.compact.Nbox)/2; %--Adjust the sub-array location of the influence function for the added zero padding

    if(any(mp.dm_ind==2)); DM2surf = padOrCropEven(DM2surf,mp.dm1.compact.NdmPad);  else; DM2surf = zeros(mp.dm1.compact.NdmPad); end 
    if(mp.flagDM2stop); DM2stop = padOrCropEven(DM2stop,mp.dm1.compact.NdmPad); else; DM2stop = ones(mp.dm1.compact.NdmPad); end
    apodRot180 = padOrCropEven( apodRot180, mp.dm1.compact.NdmPad);

    Edm1pad = padOrCropEven(Edm1,mp.dm1.compact.NdmPad); %--Pad or crop for expected sub-array indexing

    %--Propagate each actuator from DM1 through the optical system
    Gindex = 1; % initialize index counter
    for iact=mp.dm1.act_ele(:).'
        if(any(any(mp.dm1.compact.inf_datacube(:,:,iact))) && any(mp.dm1.act_ele==iact) )  %--Only compute for acutators specified for use or for influence functions that are not zeroed out
            
            %--x- and y- coordinates of the padded influence function in the full padded pupil
            x_box_AS_ind = mp.dm1.compact.xy_box_lowerLeft_AS(1,iact):mp.dm1.compact.xy_box_lowerLeft_AS(1,iact)+NboxPad1AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = mp.dm1.compact.xy_box_lowerLeft_AS(2,iact):mp.dm1.compact.xy_box_lowerLeft_AS(2,iact)+NboxPad1AS-1; % y-indices in pupil arrays for the box

            %--Propagate from DM1 to DM2, and then back to P2
            dEbox = (mirrorFac*2*pi*1j/lambda)*padOrCropEven(mp.dm1.VtoH(iact)*mp.dm1.compact.inf_datacube(:,:,iact),NboxPad1AS); %--Pad influence function at DM1 for angular spectrum propagation.
            dEbox = propcustom_PTP(dEbox.*Edm1pad(y_box_AS_ind,x_box_AS_ind),mp.P2.compact.dx*NboxPad1AS,lambda,mp.d_dm1_dm2); % forward propagate to DM2 and apply DM2 E-field
            dEP2box = propcustom_PTP(dEbox.*DM2stop(y_box_AS_ind,x_box_AS_ind).*exp(mirrorFac*2*pi*1j/lambda*DM2surf(y_box_AS_ind,x_box_AS_ind)),mp.P2.compact.dx*NboxPad1AS,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1) ); % back-propagate to DM1
            dEP2box = padOrCropEven(dEP2box,Nbox1); %--Crop down from the array size that is a power of 2 to make the MFT faster

            %--x- and y- coordinates of the UN-padded influence function in the full padded pupil
            x_box_ind = mp.dm1.compact.xy_box_lowerLeft(1,iact):mp.dm1.compact.xy_box_lowerLeft(1,iact)+Nbox1-1; % x-indices in pupil arrays for the box
            y_box_ind = mp.dm1.compact.xy_box_lowerLeft(2,iact):mp.dm1.compact.xy_box_lowerLeft(2,iact)+Nbox1-1; % y-indices in pupil arrays for the box
            x_box = mp.dm1.compact.x_pupPad(x_box_ind).'; % full pupil x-coordinates of the box 
            y_box = mp.dm1.compact.y_pupPad(y_box_ind); % full pupil y-coordinates of the box

            %--To simulate going forward to the next pupil plane (with the apodizer) most efficiently, 
            % First, back-propagate the apodizer (by rotating 180-degrees) to the previous pupil.
            % Second, negate the coordinates of the box used.
            dEP2box = apodRot180(y_box_ind,x_box_ind).*dEP2box; %--Apply 180deg-rotated SP mask.
            dEP3box = (1/1j)^2*rot90(dEP2box,2); %--Forward propagate the cropped box by rotating 180 degrees.
            x_box = rot90(-x_box,2); %--Negate to effectively rotate by 180 degrees
            y_box = rot90(-y_box,2); %--Negate to effectively rotate by 180 degrees
           
            %--Matrices for the MFT from the pupil P3 to the focal plane mask
            rect_mat_pre = (exp(-2*pi*1j*(mp.F3.compact.etas*y_box)/(lambda*mp.fl)))...
                *sqrt(mp.P2.compact.dx*mp.P2.compact.dx)*sqrt(mp.F3.compact.dxi*mp.F3.compact.deta)/(1j*lambda*mp.fl);
            rect_mat_post  = (exp(-2*pi*1j*(x_box*mp.F3.compact.xis)/(lambda*mp.fl)));

            %--MFT from pupil P3 to FPM
            EF3 = rect_mat_pre*dEP3box*rect_mat_post; % MFT to FPM
            FPM = mp.F3.compact.mask.amp.*exp(1i*2*pi/lambda*(mp.F3.n(lambda)-1)*mp.F3.t);
            EF3 = (1-FPM).*EF3; %--Apply (1-FPM) for Babinet's principle later

            %--DFT to LS ("Sub" name for Subtrahend part of the Lyot-plane E-field)
            EP4sub = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.compact.dxi,mp.F3.compact.deta,mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering);  %--Subtrahend term for the Lyot plane E-field    
            
            if(mp.useGPU);dEP2box = gather(dEP2box);end
            
            %--Full Lyot plane pupil (for Babinet)
            EP4noFPM = zeros(mp.dm1.compact.NdmPad);
            EP4noFPM(y_box_ind,x_box_ind) = dEP2box; %--Propagating the E-field from P2 to P4 without masks gives the same E-field. 
            EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.compact.Narr);
            EP4 = mp.P4.compact.croppedMask.*(EP4noFPM - EP4sub); % Babinet's principle to get E-field at Lyot plane

            % DFT to camera
            EF4 = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.F4.dxi,mp.F4.Nxi,mp.F4.deta,mp.F4.Neta,mp.centering);
            
            if(mp.useGPU);EF4 = gather(EF4);end
            
            Gzdl(:,Gindex) = EF4(mp.F4.corr.inds)/sqrt(mp.F4.compact.I00(modvar.sbpIndex));
        end
        Gindex = Gindex + 1;
    end

end    

%--DM2---------------------------------------------------------
if(whichDM==2)
    Gzdl = zeros(mp.F4.corr.Npix,mp.dm2.Nele);
    
    %--Two array sizes (at same resolution) of influence functions for MFT and angular spectrum
    Nbox2 = mp.dm2.compact.Nbox;
    NboxPad2AS = mp.dm2.compact.NboxAS; 
    mp.dm2.compact.xy_box_lowerLeft_AS = mp.dm2.compact.xy_box_lowerLeft - (NboxPad2AS-Nbox2)/2; %--Account for the padding of the influence function boxes
    
    apodRot180 = padOrCropEven( apodRot180, mp.dm2.compact.NdmPad);
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
            dEP2box = propcustom_PTP(dEbox.*Edm2(y_box_AS_ind,x_box_AS_ind),mp.P2.compact.dx*NboxPad2AS,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1) ); % back-propagate to pupil P2
            dEP2box = padOrCropEven(dEP2box,Nbox2); %--Crop down from the array size that is a power of 2 to make the MFT faster

            %--x- and y- coordinates of the UN-padded influence function in the full padded pupil
            x_box_ind = mp.dm2.compact.xy_box_lowerLeft(1,iact):mp.dm2.compact.xy_box_lowerLeft(1,iact)+Nbox2-1; % x-indices in pupil arrays for the box
            y_box_ind = mp.dm2.compact.xy_box_lowerLeft(2,iact):mp.dm2.compact.xy_box_lowerLeft(2,iact)+Nbox2-1; % y-indices in pupil arrays for the box
            x_box = mp.dm2.compact.x_pupPad(x_box_ind).'; % full pupil x-coordinates of the box 
            y_box = mp.dm2.compact.y_pupPad(y_box_ind); % full pupil y-coordinates of the box 

            %--To simulate going forward to the next pupil plane (with the apodizer) most efficiently, 
            % First, back-propagate the apodizer (by rotating 180-degrees) to the previous pupil.
            % Second, negate the coordinates of the box used.
            dEP2box = apodRot180(y_box_ind,x_box_ind).*dEP2box; %--Apply 180deg-rotated SP mask.
            dEP3box = (1/1j)^2*rot90(dEP2box,2); %--Forward propagate the cropped box by rotating 180 degrees.
            x_box = rot90(-x_box,2); %--Negate to effectively rotate by 180 degrees
            y_box = rot90(-y_box,2); %--Negate to effectively rotate by 180 degrees
            
            %--Matrices for the MFT from the pupil P3 to the focal plane mask
            rect_mat_pre = (exp(-2*pi*1j*(mp.F3.compact.etas*y_box)/(lambda*mp.fl)))...
                *sqrt(mp.P2.compact.dx*mp.P2.compact.dx)*sqrt(mp.F3.compact.dxi*mp.F3.compact.deta)/(1j*lambda*mp.fl);
            rect_mat_post  = (exp(-2*pi*1j*(x_box*mp.F3.compact.xis)/(lambda*mp.fl)));

            %--MFT from pupil P3 to FPM
            dEP2box = padOrCropEven(dEP2box,Nbox2); %--Crop back down to make the MFT faster
            EF3 = rect_mat_pre*dEP3box*rect_mat_post; % MFT to FPM
            FPM = mp.F3.compact.mask.amp.*exp(1i*2*pi/lambda*(mp.F3.n(lambda)-1)*mp.F3.t);
            EF3 = (1-FPM).*EF3; %--Apply (1-FPM) for Babinet's principle later

            % DFT to LS ("Sub" name for Subtrahend part of the Lyot-plane E-field)
            EP4sub = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.compact.dxi,mp.F3.compact.deta,mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering);  %--Subtrahend term for the Lyot plane E-field    
            
            if(mp.useGPU);dEP2box = gather(dEP2box);end
            
            EP4noFPM = zeros(mp.dm2.compact.NdmPad);
            EP4noFPM(y_box_ind,x_box_ind) = dEP2box; %--Propagating the E-field from P2 to P4 without masks gives the same E-field.
            %EP4noFPM = (1/1j)^2*rot90(EP4noFPM,2); if( strcmpi(mp.centering,'pixel') ); EP4noFPM = circshift(EP4noFPM,[1 1]); end %--Re-image to next pupil plane. (1j)^2 comes from the coefficients of the 2 skipped MFTs
            EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.compact.Narr);
            EP4 = mp.P4.compact.croppedMask.*(EP4noFPM - EP4sub); % Babinet's principle to get E-field at Lyot plane

            % DFT to camera
            EF4 = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.F4.dxi,mp.F4.Nxi,mp.F4.deta,mp.F4.Neta,mp.centering);

            if(mp.useGPU);EF4 = gather(EF4);end
            
            Gzdl(:,Gindex) = EF4(mp.F4.corr.inds)/sqrt(mp.F4.compact.I00(modvar.sbpIndex));
        end
        Gindex = Gindex + 1;
    end

end



if(mp.useGPU)
    Gzdl = gather(Gzdl);
end

end %--END OF FUNCTION


    
