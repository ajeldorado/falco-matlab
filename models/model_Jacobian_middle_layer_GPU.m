% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function jac = model_Jacobian(mp, DM, modvar)
%--Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%
% REVISION HISTORY:
% --------------
% Modified on 2018-03-06 by A.J. Riggs to be more compatible with parfor.
% Modified on 2017-11-13 by A.J. Riggs to be compatible with parfor. Added
% this function as an extra wrapper layer.
% Modified on 2017-11-09 by A.J. Riggs from model_compact.m to
%   model_Jacobian.m
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
% -vals_list = structure containing combinations of:
%     -tsi = index of the pair of sub-bandpass index and tip/tilt offset index
%     -whichDM = which DM number
%
% OUTPUTS:
% -Jac = Jacobian for the specified DM and specified T/T-wavelength pair

function jacStruct = model_Jacobian_middle_layer_GPU(mp, DM)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mirrorFac = 2; % Phase change is twice the DM surface height.
    NdmPad = DM.compact.NdmPad;

    Edm1padCube = zeros(DM.dm1.compact.NdmPad,DM.dm1.compact.NdmPad,mp.Nttlam); %--Store DM1 starting E-fields for GPU calculations
    Edm2cube = zeros(DM.dm1.compact.NdmPad,DM.dm1.compact.NdmPad,mp.Nttlam); %--Store DM2 starting E-fields for GPU calculations



for tsi=1:mp.Nttlam
    
    modvar.sbpIndex = mp.Wttlam_si(tsi);
    modvar.ttIndex = mp.Wttlam_ti(tsi);
    lambda = mp.sbp_center_vec(modvar.sbpIndex);
    
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

    if(any(DM.dm_ind==1)); DM1surf = padOrCropEven(DM.dm1.compact.surfM, NdmPad);  else; DM1surf = 0; end 
    if(any(DM.dm_ind==2)); DM2surf = padOrCropEven(DM.dm2.compact.surfM, NdmPad);  else; DM2surf = 0; end 

    if(mp.useGPU)
        pupil = gpuArray(pupil);
        Ein = gpuArray(Ein);
        if(any(DM.dm_ind==1)); DM1surf = gpuArray(DM1surf); end
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
    if(any(DM.dm_ind==1))

        DM.dm1.compact.xy_box_lowerLeft_AS = DM.dm1.compact.xy_box_lowerLeft - (DM.dm1.compact.NboxAS-DM.dm1.compact.Nbox)/2; %--Adjust the sub-array location of the influence function for the added zero padding

        if(any(DM.dm_ind==2)); DM2surf = padOrCropEven(DM2surf,DM.dm1.compact.NdmPad);  else; DM2surf = zeros(DM.dm1.compact.NdmPad); end 
        if(mp.flagDM2stop); DM2stop = padOrCropEven(DM2stop,DM.dm1.compact.NdmPad); else; DM2stop = ones(DM.dm1.compact.NdmPad); end
        apodRot180 = padOrCropEven( apodRot180, DM.dm1.compact.NdmPad);

        Edm1padCube(:,:,tsi) = padOrCropEven(Edm1,DM.dm1.compact.NdmPad); %--Pad or crop for expected sub-array indexing

    end
    %-----------------------------------------------------------

    %--DM2---------------------------------------------------------
    if(any(DM.dm_ind==2))

        DM.dm2.compact.xy_box_lowerLeft_AS = DM.dm2.compact.xy_box_lowerLeft - (DM.dm2.compact.NboxAS-DM.dm2.compact.Nbox)/2; %--Account for the padding of the influence function boxes

        apodRot180 = padOrCropEven( apodRot180, DM.dm2.compact.NdmPad);
        DM2stop = padOrCropEven(DM2stop,DM.dm2.compact.NdmPad);

        %--Propagate full field to DM2 before back-propagating in small boxes
        Edm2inc = padOrCropEven( propcustom_PTP(Edm1,DM.compact.NdmPad*mp.P2.compact.dx,lambda,mp.d_dm1_dm2), DM.dm2.compact.NdmPad); % E-field incident upon DM2
        Edm2inc = padOrCropEven(Edm2inc,DM.dm2.compact.NdmPad);
        Edm2cube(:,:,tsi) = DM2stop.*Edm2inc.*exp(mirrorFac*2*pi*1j/lambda*padOrCropEven(DM2surf,DM.dm2.compact.NdmPad)); % Initial E-field at DM2 including its own phase contribution
    end
    %-----------------------------------------------------------
   
end

    %%%%%%%%%%%%%%%
    
    %--Check which DMs are used
    if(any(DM.dm_ind==1)); gv.flagDM1 = true; else; gv.flagDM1 = false; end
    if(any(DM.dm_ind==2)); gv.flagDM2 = true; else; gv.flagDM2 = false; end
    
    
    %--"gv" structure stores all variables needed in the for loop over actuators
    gv.Edm1padCube = Edm1padCube; 
    gv.Edm2cube = Edm2cube; 
    
    gv.centering = mp.centering;
    gv.sbp_center_vec = mp.sbp_center_vec;
    gv.Nttlam = mp.Nttlam;
    gv.Wttlam_si = mp.Wttlam_si;
    gv.Wttlam_ti = mp.Wttlam_ti;
    
    gv.DM2stop = DM2stop;
    gv.DM2surf = DM2surf; 
    gv.Nbox1 = DM.dm1.compact.Nbox; 
    gv.Nbox2 = DM.dm2.compact.Nbox;     
    gv.NboxPad1AS = DM.dm1.compact.NboxAS;
    gv.NboxPad2AS = DM.dm2.compact.NboxAS; 
    gv.dm1_NdmPad = DM.dm1.compact.NdmPad;
    gv.dm2_NdmPad = DM.dm2.compact.NdmPad;
    
% %     % Move calculation of coordinate matrix, FR2shift, outside the GPU
% %     % Need a different FR2shift matrix for each different size of "Ein" array.
% %     M1 = gv.NboxPad1AS; 
% %     fx = ( -M1/2:1:(M1/2-1) )/(M1*mp.P2.compact.dx);    % frequency plane coordinates, vector
% %     [FX,FY]=meshgrid(fx); % frequency plane coordinates, matrices
% %     gv.dm1_FR2shift = fftshift(FX.^2+FY.^2); % fftshifted quadratic phase factor
% %     M2 = gv.NboxPad2AS; 
% %     fx = ( -M2/2:1:(M2/2-1) )/(M2*mp.P2.compact.dx);    % frequency plane coordinates, vector
% %     [FX,FY]=meshgrid(fx); % frequency plane coordinates, matrices
% %     gv.dm2_FR2shift = fftshift(FX.^2+FY.^2); % fftshifted quadratic phase factor
  
    gv.NactTotal1 = DM.dm1.NactTotal;
    gv.NactTotal2 = DM.dm2.NactTotal;
    gv.NactUsed1 = length(DM.dm1.act_ele);
    gv.NactUsed2 = length(DM.dm2.act_ele);
    gv.act_ele1 = DM.dm1.act_ele;
    gv.act_ele2 = DM.dm2.act_ele;
    
    gv.inf_datacube1 = DM.dm1.compact.inf_datacube;
    gv.inf_datacube2 = DM.dm2.compact.inf_datacube;    

    gv.xy_box_lowerLeft_AS_dm1 = DM.dm1.compact.xy_box_lowerLeft_AS;
    gv.xy_box_lowerLeft_AS_dm2 = DM.dm2.compact.xy_box_lowerLeft_AS;

    gv.VtoH1 = DM.dm1.VtoH;
    gv.VtoH2 = DM.dm2.VtoH;
    gv.d_dm1_dm2 = mp.d_dm1_dm2;
    gv.d_P2_dm1 = mp.d_P2_dm1;
    
    gv.dm1_xy_box_lowerLeft = DM.dm1.compact.xy_box_lowerLeft;
    gv.dm2_xy_box_lowerLeft = DM.dm2.compact.xy_box_lowerLeft;
    gv.dm1_xs = DM.dm1.compact.x_pupPad; % x 1
    gv.dm1_ys = DM.dm1.compact.y_pupPad; % y 1
    gv.dm2_xs = DM.dm2.compact.x_pupPad; % x 2
    gv.dm2_ys = DM.dm2.compact.y_pupPad; % y 2

    gv.fl = mp.fl;
    gv.F3_xis = mp.F3.compact.xis;
    gv.F3_etas = mp.F3.compact.etas; 
    gv.F3_dxi = mp.F3.compact.dxi;
    gv.F3_deta = mp.F3.compact.deta;
    
    gv.apodRot180 = apodRot180; 
    gv.F3_mask_amp = mp.F3.compact.mask.amp;
    
    gv.P2_dx = mp.P2.compact.dx;
    gv.P4_dx = mp.P4.compact.dx;
    gv.P4_Narr = mp.P4.compact.Narr;
    
% %     gv.dm1_NdmPad = DM.dm1.compact.NdmPad;
    gv.P4_croppedMask = mp.P4.compact.croppedMask;
    
    
    gv.F4_Nxi = mp.F4.compact.Nxi;
    gv.F4_Neta = mp.F4.compact.Neta;
    gv.F4_dxi = mp.F4.compact.dxi;
    gv.F4_deta = mp.F4.compact.deta;
    
    gv.F4_corr_inds = mp.F4.compact.corr.inds; 
    gv.Ncor = mp.F4.compact.Ncorr;
    gv.F4_I00 = mp.F4.compact.I00; 
    
    
    
    
    switch mp.coro 
        
        case{'LC','DMLC','APLC'} %--pupil+DMs, occulting spot FPM, and LS.
            jacStruct = model_Jacobian_LC_GPU(gv); 
            
%         case{'SPLC','FLC'} %--Optional apodizer, binary-amplitude FPM with outer diaphragm, LS
%             Jac  = model_Jacobian_SPLC(mp, DM, tsi, whichDM); 
%             
%         case{'vortex','Vortex','VC','AVC'} %--Optional apodizer, vortex FPM, LS
%             Jac  = model_Jacobian_VC(mp, DM, tsi, whichDM); 
%             
%         %case{'SPC','APP','APC'} %--Pupil-plane apodizer is only coronagraphic mask
%             %Jac  = model_Jacobian_APC(mp, DM, tsi, whichDM); 
            
        otherwise
            error('model_Jacobian_middle_layer: CASE NOT RECOGNIZED IN model_Jacobian.m');        
    end    

end
