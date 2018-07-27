% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function jac = model_Jacobian_LC(mp, DM, tsi, whichDM)
%--Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%  This function is for the Lyot coronagraph, DMLC, and APLC.
%
% REVISION HISTORY:
% --------------
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
%  gv = minimal structure of variables to send to the Jacobian calculation
%
% OUTPUTS:
% -Gttlam = Jacobian for the specified DM and specified T/T-wavelength pair
%

function jacStruc = model_Jacobian_LC_GPU(gv)

mirrorFac = 2; % Phase change is twice the DM surface height.

%--Initialize the Jacobian arrays
if(gv.flagDM1); jacStruc.G1 = double( zeros(gv.Ncor, gv.NactUsed1, gv.Nttlam) ); end % control Jacobian for DM1
if(gv.flagDM2); jacStruc.G2 = double( zeros(gv.Ncor, gv.NactUsed2 ,gv.Nttlam) ); end % control Jacobian for DM2
% if(gv.flagDM1); jacStruc.G1 = single( zeros(gv.Ncor, gv.NactUsed1, gv.Nttlam) ); end % control Jacobian for DM1
% if(gv.flagDM2); jacStruc.G2 = single( zeros(gv.Ncor, gv.NactUsed2 ,gv.Nttlam) ); end % control Jacobian for DM2




for imode = 1:gv.Nttlam
    sbpIndex = gv.Wttlam_si(imode); %--Which index in wavelength vector
    % ttIndex = gv.Wttlam_ti(imode);  %--Which index in tip/tilt offset vector  
    lambda = gv.sbp_center_vec(sbpIndex);

    %--------------------------- DM1  -------------------------------
    if(gv.flagDM1)
        Edm1pad = gv.Edm1padCube(:,:,imode);
        
        for iarray=1:gv.NactUsed1 %--Loop over DM actuators
            %--Propagate each actuator from DM1 through the optical system:
            iact = gv.act_ele1(iarray);

            %--x- and y- coordinates of the padded influence function in the full padded pupil
            x_box_AS_ind = gv.xy_box_lowerLeft_AS_dm1(1,iact):gv.xy_box_lowerLeft_AS_dm1(1,iact)+gv.NboxPad1AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = gv.xy_box_lowerLeft_AS_dm1(2,iact):gv.xy_box_lowerLeft_AS_dm1(2,iact)+gv.NboxPad1AS-1; % y-indices in pupil arrays for the box

            %--Propagate from DM1 to DM2, and then back to P2
            dEbox = (mirrorFac*2*pi*1j/lambda)*padOrCropEven(gv.VtoH1(iact)*gv.inf_datacube1(:,:,iact),gv.NboxPad1AS); %--Pad influence function at DM1 for angular spectrum propagation.
            dEbox = propcustom_PTP(dEbox.*Edm1pad(y_box_AS_ind,x_box_AS_ind),gv.P2_dx*gv.NboxPad1AS,lambda,gv.d_dm1_dm2); % forward propagate to DM2 and apply DM2 E-field
            dEP2box = propcustom_PTP(dEbox.*gv.DM2stop(y_box_AS_ind,x_box_AS_ind).*exp(mirrorFac*2*pi*1j/lambda*gv.DM2surf(y_box_AS_ind,x_box_AS_ind)),gv.P2_dx*gv.NboxPad1AS,lambda,-1*(gv.d_dm1_dm2 + gv.d_P2_dm1) ); % back-propagate to DM1
            dEP2box = padOrCropEven(dEP2box,gv.Nbox1); %--Crop down from the array size that is a power of 2 to make the MFT faster

            %--x- and y- coordinates of the UN-padded influence function in the full padded pupil
            x_box_ind =  gv.dm1_xy_box_lowerLeft(1,iact): gv.dm1_xy_box_lowerLeft(1,iact)+gv.Nbox1-1; % x-indices in pupil arrays for the box
            y_box_ind =  gv.dm1_xy_box_lowerLeft(2,iact): gv.dm1_xy_box_lowerLeft(2,iact)+gv.Nbox1-1; % y-indices in pupil arrays for the box
            x_box = gv.dm1_xs(x_box_ind).'; % full pupil x-coordinates of the box 
            y_box = gv.dm1_ys(y_box_ind); % full pupil y-coordinates of the box

            %--To simulate going forward to the next pupil plane (with the apodizer) most efficiently, 
            % First, back-propagate the apodizer (by rotating 180-degrees) to the previous pupil.
            % Second, negate the coordinates of the box used.
            dEP2box = gv.apodRot180(y_box_ind,x_box_ind).*dEP2box; %--Apply 180deg-rotated SP mask.
            dEP3box = (1/1j)^2*rot90(dEP2box,2); %--Forward propagate the cropped box by rotating 180 degrees.
            x_box = rot90(-x_box,2); %--Negate to effectively rotate by 180 degrees
            y_box = rot90(-y_box,2); %--Negate to effectively rotate by 180 degrees

            %--Matrices for the MFT from the pupil P3 to the focal plane mask
            rect_mat_pre = (exp(-2*pi*1j*(gv.F3_etas*y_box)/(lambda*gv.fl)))...
                *sqrt(gv.P2_dx*gv.P2_dx)*sqrt(gv.F3_dxi*gv.F3_deta)/(1j*lambda*gv.fl);
            rect_mat_post  = (exp(-2*pi*1j*(x_box*gv.F3_xis)/(lambda*gv.fl)));

            %--MFT from pupil P3 to FPM
            EF3 = rect_mat_pre*dEP3box*rect_mat_post; % MFT to FPM
            EF3 = (1-gv.F3_mask_amp).*EF3; %--Propagate through (1-complex FPM) for Babinet's principle

            %--DFT to LS ("Sub" name for Subtrahend part of the Lyot-plane E-field)
            EP4sub = propcustom_mft_FtoP(EF3,gv.fl,lambda,gv.F3_dxi,gv.F3_deta,gv.P4_dx,gv.P4_Narr,gv.centering);  %--Subtrahend term for the Lyot plane E-field    

            %--Full Lyot plane pupil (for Babinet)
            EP4noFPM = zeros(gv.dm1_NdmPad);
            EP4noFPM(y_box_ind,x_box_ind) = dEP2box; %--Propagating the E-field from P2 to P4 without masks gives the same E-field. 
            EP4noFPM = padOrCropEven(EP4noFPM,gv.P4_Narr);
            EP4 = gv.P4_croppedMask.*(EP4noFPM - EP4sub); % Babinet's principle to get E-field at Lyot plane

            % DFT to camera
            EF4 = propcustom_mft_PtoF(EP4,gv.fl,lambda,gv.P4_dx,gv.F4_dxi,gv.F4_Nxi,gv.F4_deta,gv.F4_Neta,gv.centering);

            jacStruc.G1(:,iarray,imode) = EF4(gv.F4_corr_inds)/sqrt(gv.F4_I00(sbpIndex));
            
        end

    end    



    %--------------------------- DM2  -------------------------------
    if(gv.flagDM2)
        Edm2 = gv.Edm2cube(:,:,imode);

        for iarray=1:gv.NactUsed2 %--Loop over DM actuators
            %--Propagate each actuator from DM2 through the rest of the optical system:
            iact = gv.act_ele2(iarray);

            %--x- and y- coordinates of the padded influence function in the full padded pupil
            x_box_AS_ind = gv.xy_box_lowerLeft_AS_dm2(1,iact):gv.xy_box_lowerLeft_AS_dm2(1,iact)+gv.NboxPad2AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = gv.xy_box_lowerLeft_AS_dm2(2,iact):gv.xy_box_lowerLeft_AS_dm2(2,iact)+gv.NboxPad2AS-1; % y-indices in pupil arrays for the box

            dEbox = gv.VtoH2(iact)*(mirrorFac*2*pi*1j/lambda)*padOrCropEven(gv.inf_datacube2(:,:,iact),gv.NboxPad2AS); %--the padded influence function at DM2
            dEP2box = propcustom_PTP(dEbox.*Edm2(y_box_AS_ind,x_box_AS_ind),gv.P2_dx*gv.NboxPad2AS,lambda,-1*(gv.d_dm1_dm2 + gv.d_P2_dm1) ); % back-propagate to pupil P2
            dEP2box = padOrCropEven(dEP2box,gv.Nbox2); %--Crop down from the array size that is a power of 2 to make the MFT faster

            %--x- and y- coordinates of the UN-padded influence function in the full padded pupil
            x_box_ind =  gv.dm2_xy_box_lowerLeft(1,iact): gv.dm2_xy_box_lowerLeft(1,iact)+gv.Nbox2-1; % x-indices in pupil arrays for the box
            y_box_ind =  gv.dm2_xy_box_lowerLeft(2,iact): gv.dm2_xy_box_lowerLeft(2,iact)+gv.Nbox2-1; % y-indices in pupil arrays for the box
            x_box = gv.dm2_xs(x_box_ind).'; % full pupil x-coordinates of the box 
            y_box = gv.dm2_ys(y_box_ind); % full pupil y-coordinates of the box 

            %--To simulate going forward to the next pupil plane (with the apodizer) most efficiently, 
            % First, back-propagate the apodizer (by rotating 180-degrees) to the previous pupil.
            % Second, negate the coordinates of the box used.
            dEP2box = gv.apodRot180(y_box_ind,x_box_ind).*dEP2box; %--Apply 180deg-rotated SP mask.
            dEP3box = (1/1j)^2*rot90(dEP2box,2); %--Forward propagate the cropped box by rotating 180 degrees.
            x_box = rot90(-x_box,2); %--Negate to effectively rotate by 180 degrees
            y_box = rot90(-y_box,2); %--Negate to effectively rotate by 180 degrees

            %--Matrices for the MFT from the pupil P3 to the focal plane mask
            rect_mat_pre = (exp(-2*pi*1j*(gv.F3_etas*y_box)/(lambda*gv.fl)))...
                *sqrt(gv.P2_dx*gv.P2_dx)*sqrt(gv.F3_dxi*gv.F3_deta)/(1j*lambda*gv.fl);
            rect_mat_post  = (exp(-2*pi*1j*(x_box*gv.F3_xis)/(lambda*gv.fl)));

            %--MFT from pupil P3 to FPM
            dEP2box = padOrCropEven(dEP2box,gv.Nbox2); %--Crop back down to make the MFT faster
            EF3 = rect_mat_pre*dEP3box*rect_mat_post; % MFT to FPM
            EF3 = (1-gv.F3_mask_amp).*EF3; %--Propagate through ( 1 - (complex FPM) ) for Babinet's principle

            % DFT to LS ("Sub" name for Subtrahend part of the Lyot-plane E-field)
            EP4sub = propcustom_mft_FtoP(EF3,gv.fl,lambda,gv.F3_dxi,gv.F3_deta,gv.P4_dx,gv.P4_Narr,gv.centering);  %--Subtrahend term for the Lyot plane E-field    

            EP4noFPM = zeros(gv.dm2_NdmPad);
            EP4noFPM(y_box_ind,x_box_ind) = dEP2box; %--Propagating the E-field from P2 to P4 without masks gives the same E-field.
            EP4noFPM = padOrCropEven(EP4noFPM,gv.P4_Narr);
            EP4 = gv.P4_croppedMask.*(EP4noFPM - EP4sub); % Babinet's principle to get E-field at Lyot plane

            % DFT to camera
            EF4 = propcustom_mft_PtoF(EP4,gv.fl,lambda,gv.P4_dx,gv.F4_dxi,gv.F4_Nxi,gv.F4_deta,gv.F4_Neta,gv.centering);

            jacStruc.G2(:,iarray,imode) = EF4(gv.F4_corr_inds)/sqrt(gv.F4_I00(sbpIndex));
            
        end

    end


end


end %--END OF FUNCTION


    
