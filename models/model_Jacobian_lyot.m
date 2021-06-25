% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%
%  For lyot-type coronagraphs with or without an outer stop on the FPM, 
%  with or without an apodizer, and with an opaque or complex-valued
%  occulting spot.
%
% INPUTS
% ------
% mp : structure of model parameters
% iMode : index of the pair of sub-bandpass index and Zernike mode index
% whichDM : which DM number
%
% OUTPUTS
% ------
% Gmode : Jacobian for the specified Zernike mode, DM, star, and subband.

function Gmode = model_Jacobian_lyot(mp, iMode, whichDM)

modvar.sbpIndex = mp.jac.sbp_inds(iMode);
modvar.zernIndex = mp.jac.zern_inds(iMode);
modvar.starIndex = mp.jac.star_inds(iMode);
lambda = mp.sbp_centers(modvar.sbpIndex); 
NdmPad = mp.compact.NdmPad;
surfIntoPhase = 2;
scaleFac = 1; % Default is that F3 focal plane sampling does not vary with wavelength

if mp.flagRotation
    NrelayFactor = 1;
else
    NrelayFactor = 0; % zero out the number of relays
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input E-field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Include the star position and weight in the starting wavefront
iStar = modvar.starIndex;
xiOffset = mp.compact.star.xiOffsetVec(iStar);
etaOffset = mp.compact.star.etaOffsetVec(iStar);
starWeight = mp.compact.star.weights(iStar);
TTphase = (-1)*(2*pi*(xiOffset*mp.P2.compact.XsDL + etaOffset*mp.P2.compact.YsDL));
Ett = exp(1j*TTphase*mp.lambda0/lambda);
Ein = sqrt(starWeight) * Ett .* mp.P1.compact.E(:, :, modvar.sbpIndex);

%--Apply a Zernike (in amplitude) at input pupil
if modvar.zernIndex ~= 1
    indsZnoll = modvar.zernIndex; %--Just send in 1 Zernike mode
    zernMat = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam, mp.centering, indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    zernMat = pad_crop(zernMat, mp.P1.compact.Narr);
    Ein = Ein .* zernMat * (2*pi*1j/lambda) * mp.jac.Zcoef(mp.jac.zerns == modvar.zernIndex);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pupil = pad_crop(mp.P1.compact.mask, NdmPad);
Ein = pad_crop(Ein, NdmPad);
if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
end

%--Re-image the apodizer from pupil P3 back to pupil P2.
if mp.flagApod
    apodReimaged = pad_crop(mp.P3.compact.mask, NdmPad);
    apodReimaged = propcustom_relay(apodReimaged, NrelayFactor*mp.Nrelay2to3, mp.centering);
else
    apodReimaged = ones(NdmPad); 
end

if mp.flagDM1stop; DM1stop = pad_crop(mp.dm1.compact.mask, NdmPad); else; DM1stop = ones(NdmPad); end
if mp.flagDM2stop; DM2stop = pad_crop(mp.dm2.compact.mask, NdmPad); else; DM2stop = ones(NdmPad); end

if any(mp.dm_ind == 1); DM1surf = pad_crop(mp.dm1.compact.surfM, NdmPad);  else; DM1surf = 0; end 
if any(mp.dm_ind == 2); DM2surf = pad_crop(mp.dm2.compact.surfM, NdmPad);  else; DM2surf = 0; end 
if mp.useGPU
    if any(mp.dm_ind == 1); DM1surf = gpuArray(DM1surf); end
    if any(mp.dm_ind == 2); DM2surf = gpuArray(DM2surf); end
end

% Define the FPM
switch upper(mp.coro)
    
    case{'LC', 'APLC', 'FLC', 'SPLC'}
        fpm = mp.F3.compact.mask;
        transOuterFPM = 1;
        
    case{'HLC'}
        switch mp.layout
            case{'fourier'}
                %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
                t_Ti_base = 0;
                t_Ni_vec = 0;
                t_PMGI_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
                pol = 2;
                [transOuterFPM, ~] = falco_thin_film_material_def(lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, lambda*mp.FPM.d0fac, pol);
                fpm = squeeze(mp.FPMcube(:, :, modvar.sbpIndex)); %--Complex transmission of the FPM. Calculated in model_Jacobian.m
            case{'fpm_scale', 'proper', 'roman_phasec_proper', 'wfirst_phaseb_proper'}
                fpm = squeeze(mp.compact.FPMcube(:, :, modvar.sbpIndex)); %--Complex transmission of the FPM.
                transOuterFPM = fpm(1, 1); %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
                scaleFac = lambda/mp.lambda0; % Focal plane sampling varies with wavelength
            otherwise
                error('Invalid combination of mp.layout and mp.coro')
        end
        
    otherwise
        error('Value of mp.coro not recognized.');
end

xisF3 = scaleFac * mp.F3.compact.xis;
etasF3 = scaleFac * mp.F3.compact.etas;

%--For including DM surface errors (quilting, scalloping, etc.)
if mp.flagDMwfe
    if any(mp.dm_ind == 1); Edm1WFE = exp(2*pi*1j/lambda.*pad_crop(mp.dm1.compact.wfe, NdmPad, 'extrapval', 0)); else; Edm1WFE = ones(NdmPad); end
    if any(mp.dm_ind == 2); Edm2WFE = exp(2*pi*1j/lambda.*pad_crop(mp.dm2.compact.wfe, NdmPad, 'extrapval', 0)); else; Edm2WFE = ones(NdmPad); end
else
    Edm1WFE = ones(NdmPad);
    Edm2WFE = ones(NdmPad);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the unocculted peak E-field and coronagraphic E-field
if mp.jac.minimizeNI
    modvar.whichSource = 'star';
    Eocculted = model_compact(mp, modvar);
    Eunocculted = model_compact(mp, modvar, 'nofpm');
    [~, indPeak] = max(abs(Eunocculted(:)));
    Epeak = Eunocculted(indPeak);
end

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil .* Ein; %--E-field at pupil plane P1
EP2 = propcustom_relay(EP1, NrelayFactor*mp.Nrelay1to2, mp.centering); 

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
Edm1 = propcustom_PTP(EP2, mp.P2.compact.dx*NdmPad, lambda, mp.d_P2_dm1);
Edm1 = Edm1WFE .* DM1stop .* exp(surfIntoPhase*2*pi*1j*DM1surf/lambda) .* Edm1;

%--DM1---------------------------------------------------------
if whichDM == 1 
    Gmode = zeros(mp.Fend.corr.Npix, mp.dm1.Nele);
    
    %--Two array sizes (at same resolution) of influence functions for MFT and angular spectrum
    NboxPad1AS = mp.dm1.compact.NboxAS;
    mp.dm1.compact.xy_box_lowerLeft_AS = mp.dm1.compact.xy_box_lowerLeft - (mp.dm1.compact.NboxAS-mp.dm1.compact.Nbox)/2; %--Adjust the sub-array location of the influence function for the added zero padding

    if any(mp.dm_ind == 2)
        DM2surf = pad_crop(DM2surf, mp.dm1.compact.NdmPad);
    else
        DM2surf = zeros(mp.dm1.compact.NdmPad);
    end
    
    if mp.flagDM2stop
        DM2stop = pad_crop(DM2stop, mp.dm1.compact.NdmPad);
    else
        DM2stop = ones(mp.dm1.compact.NdmPad);
    end
    
    apodReimaged = pad_crop(apodReimaged, mp.dm1.compact.NdmPad);

    Edm1pad = pad_crop(Edm1, mp.dm1.compact.NdmPad); %--Pad or crop for expected sub-array indexing
    Edm2WFEpad = pad_crop(Edm2WFE, mp.dm1.compact.NdmPad); %--Pad or crop for expected sub-array indexing
    
    %--Propagate each actuator from DM1 through the optical system
    Gindex = 1; % initialize index counter
    for iact = mp.dm1.act_ele(:).'  %--MUST BE A COLUMN VECTOR `
        if any(any(mp.dm1.compact.inf_datacube(:, :, iact)))
            %--x- and y- coordinates' indices of the padded influence function in the full padded pupil
            x_box_AS_ind = mp.dm1.compact.xy_box_lowerLeft_AS(1, iact):mp.dm1.compact.xy_box_lowerLeft_AS(1, iact)+NboxPad1AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = mp.dm1.compact.xy_box_lowerLeft_AS(2, iact):mp.dm1.compact.xy_box_lowerLeft_AS(2, iact)+NboxPad1AS-1; % y-indices in pupil arrays for the box
            %--x- and y- coordinates of the padded influence function in the full padded pupil
            x_box = mp.dm1.compact.x_pupPad(x_box_AS_ind).'; % full pupil x-coordinates of the box 
            y_box = mp.dm1.compact.y_pupPad(y_box_AS_ind); % full pupil y-coordinates of the box
            
            %--Propagate from DM1 to DM2, and then back to P2
            dEbox = (surfIntoPhase*2*pi*1j/lambda)*pad_crop(mp.dm1.VtoH(iact)*mp.dm1.compact.inf_datacube(:, :, iact), NboxPad1AS); %--Pad influence function at DM1 for angular spectrum propagation.
            dEbox = propcustom_PTP_inf_func(dEbox.*Edm1pad(y_box_AS_ind, x_box_AS_ind), mp.P2.compact.dx*NboxPad1AS, lambda, mp.d_dm1_dm2, mp.dm1.dm_spacing, mp.propMethodPTP); % forward propagate to DM2 and apply DM2 E-field
            dEP2box = propcustom_PTP_inf_func(dEbox.*Edm2WFEpad(y_box_AS_ind, x_box_AS_ind).*DM2stop(y_box_AS_ind, x_box_AS_ind).*exp(surfIntoPhase*2*pi*1j/lambda*DM2surf(y_box_AS_ind, x_box_AS_ind)), mp.P2.compact.dx*NboxPad1AS, lambda, -1*(mp.d_dm1_dm2 + mp.d_P2_dm1), mp.dm1.dm_spacing, mp.propMethodPTP ); % back-propagate to DM1
            
            %--To simulate going forward to the next pupil plane (with the apodizer) most efficiently, 
            % First, back-propagate the apodizer (by rotating 180-degrees) to the previous pupil.
            % Second, negate the coordinates of the box used.
            dEP2box = apodReimaged(y_box_AS_ind, x_box_AS_ind) .* dEP2box;
            dEP3box = rot90(dEP2box, NrelayFactor*2*mp.Nrelay2to3);
            x_box = (-1)^(NrelayFactor*mp.Nrelay2to3) * rot90(x_box, NrelayFactor*2*mp.Nrelay2to3); %--Negate and rotate coordinates to effectively rotate by 180 degrees. No change if 360 degree rotation.
            y_box = (-1)^(NrelayFactor*mp.Nrelay2to3) * rot90(y_box, NrelayFactor*2*mp.Nrelay2to3); %--Negate and rotate coordinates to effectively rotate by 180 degrees. No change if 360 degree rotation.

            %--Matrices for the MFT from the pupil P3 to the focal plane mask
            rect_mat_pre = (exp(-2*pi*1j*(etasF3*y_box)/(lambda*mp.fl)))...
                * sqrt(mp.P2.compact.dx*mp.P2.compact.dx) * ...
                scaleFac * sqrt(mp.F3.compact.dxi*mp.F3.compact.deta) / (lambda*mp.fl);
            rect_mat_post  = (exp(-2*pi*1j*(x_box*xisF3)/(lambda*mp.fl)));

            %--MFT from pupil P3 to FPM
            EF3inc = rect_mat_pre * dEP3box * rect_mat_post; % MFT to FPM
            
            switch upper(mp.coro)
                
                case{'LC', 'APLC', 'HLC'}
                    EF3 = (transOuterFPM - fpm) .* EF3inc; %--Propagate through (1-complex FPM) for Babinet's principle

                    %--DFT to LS ("Sub" name for Subtrahend part of the Lyot-plane E-field)
                    EP4subtrahend = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.compact.dxi, ...
                        scaleFac*mp.F3.compact.deta, mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering);
                    EP4subtrahend = propcustom_relay(EP4subtrahend, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);

                    %--Full Lyot plane pupil (for Babinet)
                    EP4noFPM = zeros(mp.dm1.compact.NdmPad);
                    if mp.useGPU; EP4noFPM = gpuArray(EP4noFPM);end
                    EP4noFPM(y_box_AS_ind, x_box_AS_ind) = dEP2box; %--Propagating the E-field from P2 to P4 without masks gives the same E-field. 
                    EP4noFPM = propcustom_relay(EP4noFPM, NrelayFactor*(mp.Nrelay2to3+mp.Nrelay3to4), mp.centering); %--Get the correct orientation 
                    EP4noFPM = pad_crop(EP4noFPM, mp.P4.compact.Narr);
                    EP4 = transOuterFPM*EP4noFPM - EP4subtrahend; % Babinet's principle to get E-field at Lyot plane
                    
                case{'FLC', 'SPLC'}
                    EF3 = mp.F3.compact.mask .* EF3inc;
                    EP4 = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.compact.dxi, scaleFac*mp.F3.compact.deta, mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering); % MFT to Lyot Plane
                    EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);
                
                otherwise
                    error('Value of mp.coro not recognized.')
                    
            end

            EP4 = mp.P4.compact.croppedMask .* EP4; % Apply Lyot stop
            
            %--MFT to final focal plane
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.NrelayFend, mp.centering);
            EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, mp.Fend.dxi, ...
                mp.Fend.Nxi, mp.Fend.deta, mp.Fend.Neta, mp.centering);
            if(mp.useGPU); EFend = gather(EFend); end
            
            Gmode(:, Gindex) = EFend(mp.Fend.corr.maskBool) / sqrt(mp.Fend.compact.I00(modvar.sbpIndex));
        end
        
        Gindex = Gindex + 1;
    end
    
    if mp.jac.minimizeNI
       JacOfPeak = model_Jacobian_no_FPM(mp, iMode, whichDM); 
       Gmode = Gmode/Epeak - Eocculted(mp.Fend.corr.maskBool) / (Epeak*Epeak) .* repmat(JacOfPeak, [mp.Fend.corr.Npix, 1]);
    end
    
    Gmode = mp.dm1.weight * Gmode;
end    

%--DM2---------------------------------------------------------
if whichDM == 2
    Gmode = zeros(mp.Fend.corr.Npix, mp.dm2.Nele);
    
    %--Two array sizes (at same resolution) of influence functions for MFT and angular spectrum
    NboxPad2AS = mp.dm2.compact.NboxAS; 
    mp.dm2.compact.xy_box_lowerLeft_AS = mp.dm2.compact.xy_box_lowerLeft - (NboxPad2AS-mp.dm2.compact.Nbox)/2; %--Account for the padding of the influence function boxes
    
    apodReimaged = pad_crop(apodReimaged, mp.dm2.compact.NdmPad);
    DM2stopPad = pad_crop(DM2stop, mp.dm2.compact.NdmPad);
    Edm2WFEpad = pad_crop(Edm2WFE, mp.dm2.compact.NdmPad);
    
    %--Propagate full field to DM2 before back-propagating in small boxes
    Edm2inc = propcustom_PTP(Edm1, mp.compact.NdmPad*mp.P2.compact.dx, lambda, mp.d_dm1_dm2); % E-field incident upon DM2
    Edm2inc = pad_crop(Edm2inc, mp.dm2.compact.NdmPad);
    Edm2 = Edm2WFEpad.*DM2stopPad.*Edm2inc.*exp(surfIntoPhase*2*pi*1j/lambda*pad_crop(DM2surf, mp.dm2.compact.NdmPad)); % Initial E-field at DM2 including its own phase contribution
    
    %--Propagate each actuator from DM2 through the rest of the optical system
    Gindex = 1; % Initialize index counter
    for iact = mp.dm2.act_ele(:).'  %--Only compute for acutators specified %--MUST BE A COLUMN VECTOR
        if any(any(mp.dm2.compact.inf_datacube(:, :, iact)))
            %--x- and y- coordinates and their indices of the padded influence function in the full padded pupil
            x_box_AS_ind = mp.dm2.compact.xy_box_lowerLeft_AS(1, iact):mp.dm2.compact.xy_box_lowerLeft_AS(1, iact)+NboxPad2AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = mp.dm2.compact.xy_box_lowerLeft_AS(2, iact):mp.dm2.compact.xy_box_lowerLeft_AS(2, iact)+NboxPad2AS-1; % y-indices in pupil arrays for the box
            x_box = mp.dm2.compact.x_pupPad(x_box_AS_ind).'; % full pupil x-coordinates of the box 
            y_box = mp.dm2.compact.y_pupPad(y_box_AS_ind); % full pupil y-coordinates of the box 
            
            dEbox = mp.dm2.VtoH(iact)*(surfIntoPhase*2*pi*1j/lambda)*pad_crop(mp.dm2.compact.inf_datacube(:, :, iact), NboxPad2AS); %--the padded influence function at DM2
            dEP2box = propcustom_PTP_inf_func(dEbox.*Edm2(y_box_AS_ind, x_box_AS_ind), mp.P2.compact.dx*NboxPad2AS, lambda, -1*(mp.d_dm1_dm2 + mp.d_P2_dm1), mp.dm2.dm_spacing, mp.propMethodPTP); % back-propagate to pupil P2

            %--To simulate going forward to the next pupil plane (with the apodizer) most efficiently, 
            % First, back-propagate the apodizer (by rotating 180-degrees) to the previous pupil.
            % Second, negate the coordinates of the box used.
            dEP2box = apodReimaged(y_box_AS_ind, x_box_AS_ind).*dEP2box;
            dEP3box = rot90(dEP2box, NrelayFactor*2*mp.Nrelay2to3);
            x_box = (-1)^(NrelayFactor*mp.Nrelay2to3) * rot90(x_box, NrelayFactor*2*mp.Nrelay2to3); %--Negate and rotate coordinates to effectively rotate by 180 degrees. No change if 360 degree rotation.
            y_box = (-1)^(NrelayFactor*mp.Nrelay2to3) * rot90(y_box, NrelayFactor*2*mp.Nrelay2to3); %--Negate and rotate coordinates to effectively rotate by 180 degrees. No change if 360 degree rotation.

            %--Matrices for the MFT from the pupil P3 to the focal plane mask
            rect_mat_pre = (exp(-2*pi*1j*(etasF3*y_box)/(lambda*mp.fl)))...
                * sqrt(mp.P2.compact.dx*mp.P2.compact.dx) * ...
                scaleFac*sqrt(mp.F3.compact.dxi*mp.F3.compact.deta) / (lambda*mp.fl);
            rect_mat_post = (exp(-2*pi*1j*(x_box*xisF3)/(lambda*mp.fl)));

            %--MFT from pupil P3 to FPM
            EF3inc = rect_mat_pre*dEP3box * rect_mat_post; % MFT to FPM
            
            switch upper(mp.coro)
                
                case{'LC', 'APLC', 'HLC'}
                    EF3 = (transOuterFPM-fpm) .* EF3inc; %--Propagate through (1 - (complex FPM)) for Babinet's principle

                    % DFT to LS
                    EP4subtrahend = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.compact.dxi, ...
                        scaleFac*mp.F3.compact.deta, mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering);
                    EP4subtrahend = propcustom_relay(EP4subtrahend, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering); %--Get the correct orientation

                    EP4noFPM = zeros(mp.dm2.compact.NdmPad);
                    if mp.useGPU; EP4noFPM = gpuArray(EP4noFPM);end
                    EP4noFPM(y_box_AS_ind, x_box_AS_ind) = dEP2box; %--Propagating the E-field from P2 to P4 without masks gives the same E-field.
                    EP4noFPM = propcustom_relay(EP4noFPM, NrelayFactor*(mp.Nrelay2to3+mp.Nrelay3to4), mp.centering); %--Get the correct orientation 
                    EP4noFPM = pad_crop(EP4noFPM, mp.P4.compact.Narr);
                    EP4 = transOuterFPM*EP4noFPM - EP4subtrahend; % Babinet's principle to get E-field at Lyot plane
            
                case{'FLC', 'SPLC'}
                    EF3 = mp.F3.compact.mask .* EF3inc; % Apply FPM
                    EP4 = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.compact.dxi, ...
                        scaleFac*mp.F3.compact.deta, mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering); % MFT to Lyot Plane
                    EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering); % Add more re-imaging relays between pupils P3 and P4 if necessary
                
                otherwise
                    error('Value of mp.coro not recognized.')
                    
            end
                    
            EP4 = mp.P4.compact.croppedMask .* EP4;
            
            %--MFT to final focal plane
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.NrelayFend, mp.centering);
            EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, mp.Fend.dxi, mp.Fend.Nxi, mp.Fend.deta, mp.Fend.Neta, mp.centering);
            if mp.useGPU; EFend = gather(EFend); end

            Gmode(:, Gindex) = EFend(mp.Fend.corr.maskBool) / sqrt(mp.Fend.compact.I00(modvar.sbpIndex));
        end
        Gindex = Gindex + 1;
    end
    
    if mp.jac.minimizeNI
       JacOfPeak = model_Jacobian_no_FPM(mp, iMode, whichDM); 
       Gmode = Gmode/Epeak - Eocculted(mp.Fend.corr.maskBool) / ...
           (Epeak*Epeak) .* repmat(JacOfPeak, [mp.Fend.corr.Npix, 1]);
    end
    
    Gmode = mp.dm2.weight * Gmode;
end

%--DM8--------------------------------------------------------- 
if whichDM == 8
    Gmode = zeros(mp.Fend.corr.Npix, mp.dm8.Nele);
    Nbox8 = mp.dm8.compact.Nbox;
    
    stepFac = 1; %--Adjust the step size in the Jacobian, then divide back out. Used for helping counteract effect of discretization.
    
    %--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
    Edm2 = Edm2WFE.*DM2stop.*exp(surfIntoPhase*2*pi*1j*DM2surf/lambda).*propcustom_PTP(Edm1, mp.P2.compact.dx*NdmPad, lambda, mp.d_dm1_dm2); % Pre-compute the initial DM2 E-field
    
    %--Back-propagate to pupil P2
    EP2eff = propcustom_PTP(Edm2, mp.P2.compact.dx*NdmPad, lambda, -1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); 
    
    %--Rotate 180 degrees mp.Nrelay2to3 times to go from pupil P2 to P3
    EP3 = propcustom_relay(EP2eff, NrelayFactor*mp.Nrelay2to3, mp.centering);

    %--Apply apodizer mask.
    if mp.flagApod
        EP3 = mp.P3.compact.mask.*pad_crop(EP3, mp.P1.compact.Narr); 
    end
    
    %--MFT from pupil P3 to FPM (at focus F3)
    EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.compact.dx, scaleFac*mp.F3.compact.dxi, ...
        mp.F3.compact.Nxi, scaleFac*mp.F3.compact.deta, mp.F3.compact.Neta, mp.centering);
    EF3inc = pad_crop(EF3inc, mp.dm8.compact.NdmPad);
    %--Coordinates for metal thickness and dielectric thickness
    DM9transIndAll = falco_discretize_FPM_surf(mp.dm9.surf, mp.t_diel_nm_vec, mp.dt_diel_nm); %--All of the mask
    
    %--Propagate each actuator from DM9 through the rest of the optical system
    Gindex = 1; % initialize index counter
    for iact = mp.dm8.act_ele(:).' %--MUST BE A COLUMN VECTOR 
         if any(any(mp.dm8.compact.inf_datacube(:, :, iact)))
            %--xi- and eta- coordinates in the full FPM portion of the focal plane
            xi_box_ind = mp.dm8.compact.xy_box_lowerLeft(1, iact):...
                mp.dm8.compact.xy_box_lowerLeft(1, iact)+mp.dm8.compact.Nbox-1;
            eta_box_ind = mp.dm8.compact.xy_box_lowerLeft(2, iact):...
                mp.dm8.compact.xy_box_lowerLeft(2, iact)+mp.dm8.compact.Nbox-1;
            xi_box = mp.dm8.compact.x_pupPad(xi_box_ind).'; % full image xi-coordinates of the box 
            eta_box = mp.dm8.compact.y_pupPad(eta_box_ind); % full image eta-coordinates of the box 

            %--Obtain values for the "poked" FPM's complex transmission (only in the sub-array where poked)
            Nxi = Nbox8;
            Neta = Nbox8;
            
            DM8surfCropNew = stepFac*mp.dm8.VtoH(iact).* mp.dm8.compact.inf_datacube(:, :, iact) + ...
                mp.dm8.surf(eta_box_ind, xi_box_ind);
            DM8transInd = falco_discretize_FPM_surf(DM8surfCropNew, mp.t_metal_nm_vec, mp.dt_metal_nm);
            DM9transInd = DM9transIndAll(eta_box_ind, xi_box_ind); %--Cropped region of the FPM.

            %--Look up table to compute complex transmission coefficient of the FPM at each pixel
            FPMpoked = zeros(Neta, Nxi); %--Initialize output array of FPM's complex transmission    
            for ix = 1:Nxi
                for iy = 1:Neta
                    ind_metal = DM8transInd(iy, ix);
                    ind_diel  = DM9transInd(iy, ix);
                    FPMpoked(iy, ix) = mp.complexTransCompact(ind_diel, ind_metal, modvar.sbpIndex);
                end
            end            
  
            dEF3box = ((transOuterFPM-FPMpoked) - (transOuterFPM-fpm(eta_box_ind, xi_box_ind))) .* ...
                EF3inc(eta_box_ind, xi_box_ind); % Delta field (in a small region) at the FPM

            %--Matrices for the MFT from the FPM stamp to the Lyot stop
            rect_mat_pre = (exp(-2*pi*1j*(mp.P4.compact.ys*eta_box)/(lambda*mp.fl)))...
                *sqrt(mp.P4.compact.dx*mp.P4.compact.dx)*scaleFac*sqrt(mp.F3.compact.dxi*mp.F3.compact.deta)/(lambda*mp.fl);
            rect_mat_post  = (exp(-2*pi*1j*(xi_box*mp.P4.compact.xs)/(lambda*mp.fl)));

            %--DFT from FPM to Lyot stop (Nominal term transOuterFPM*EP4noFPM subtracts out to 0 since it ignores the FPM change).
            EP4 = 0 - rect_mat_pre*dEF3box*rect_mat_post; % MFT from FPM (F3) to Lyot stop plane (P4)
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering); %--Get the correct orientation 
            EP4 = mp.P4.compact.croppedMask .* EP4; %--Apply Lyot stop

            %--MFT to final focal plane
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.NrelayFend, mp.centering);
            EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, mp.Fend.dxi, ...
                mp.Fend.Nxi, mp.Fend.deta, mp.Fend.Neta, mp.centering);
            if mp.useGPU; EFend = gather(EFend); end
            
            Gmode(:, Gindex) = (1/stepFac)*mp.dm8.weight*EFend(mp.Fend.corr.maskBool) / ...
                sqrt(mp.Fend.compact.I00(modvar.sbpIndex));
            Gindex = Gindex + 1;
        end
    end
end

%--DM9--------------------------------------------------------- 
if whichDM == 9
    Gmode = zeros(mp.Fend.corr.Npix, mp.dm9.Nele);
    Nbox9 = mp.dm9.compact.Nbox;
    
    % Adjust the step size in the Jacobian, then divide back out later.
    % Used for helping counteract effect of discretization.
    if isfield(mp.dm9, 'stepFac') == false
        stepFac = 20;
    else
        stepFac = mp.dm9.stepFac;
    end

    %--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
    Edm2 = Edm2WFE .* DM2stop .* exp(surfIntoPhase*2*pi*1j*DM2surf/lambda) .* ...
        propcustom_PTP(Edm1, mp.P2.compact.dx*NdmPad, lambda, mp.d_dm1_dm2);
    
    %--Back-propagate to pupil P2
    EP2eff = propcustom_PTP(Edm2, mp.P2.compact.dx*NdmPad, lambda, -1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); 
    
    %--Rotate 180 degrees mp.Nrelay2to3 times to go from pupil P2 to P3
    EP3 = propcustom_relay(EP2eff, NrelayFactor*mp.Nrelay2to3, mp.centering);

    %--Apply apodizer mask.
    if mp.flagApod
        EP3 = mp.P3.compact.mask.*pad_crop(EP3, mp.P1.compact.Narr); 
    end
    
    %--MFT from pupil P3 to FPM (at focus F3)
    EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.compact.dx, scaleFac*mp.F3.compact.dxi, ...
        mp.F3.compact.Nxi, scaleFac*mp.F3.compact.deta, mp.F3.compact.Neta, mp.centering);
    EF3inc = pad_crop(EF3inc, mp.dm9.compact.NdmPad);
    
    %--Coordinates for metal thickness and dielectric thickness
    DM8transIndAll = falco_discretize_FPM_surf(mp.dm8.surf, mp.t_metal_nm_vec, mp.dt_metal_nm); %--All of the mask
    
    %--Propagate each actuator from DM9 through the rest of the optical system
    Gindex = 1; % initialize index counter
    for iact = mp.dm9.act_ele(:).' %--MUST BE A COLUMN VECTOR
        if any(any(mp.dm9.compact.inf_datacube(:, :, iact)))
            %--xi- and eta- coordinates in the full FPM portion of the focal plane
            xi_box_ind = mp.dm9.compact.xy_box_lowerLeft(1, iact):...
                mp.dm9.compact.xy_box_lowerLeft(1, iact)+mp.dm9.compact.Nbox-1; % xi-indices in image arrays for the box
            eta_box_ind = mp.dm9.compact.xy_box_lowerLeft(2, iact):...
                mp.dm9.compact.xy_box_lowerLeft(2, iact)+mp.dm9.compact.Nbox-1; % eta-indices in image arrays for the box
            xi_box = mp.dm9.compact.x_pupPad(xi_box_ind).'; % full image xi-coordinates of the box 
            eta_box = mp.dm9.compact.y_pupPad(eta_box_ind); % full image eta-coordinates of the box 

            %--Obtain values for the "poked" FPM's complex transmission (only in the sub-array where poked)
            Nxi = Nbox9;
            Neta = Nbox9;
            DM9surfCropNew = stepFac*mp.dm9.VtoH(iact).*mp.dm9.compact.inf_datacube(:, :, iact) + ...
                mp.dm9.surf(eta_box_ind, xi_box_ind); % New DM9 surface profile in the poked region (meters)
            DM9transInd = falco_discretize_FPM_surf(DM9surfCropNew, mp.t_diel_nm_vec, mp.dt_diel_nm);
            DM8transInd = DM8transIndAll(eta_box_ind, xi_box_ind); %--Cropped region of the FPM.

            %--Look up table to compute complex transmission coefficient of the FPM at each pixel
            FPMpoked = zeros(Neta, Nxi); %--Initialize output array of FPM's complex transmission    
            for ix = 1:Nxi
                for iy = 1:Neta
                    ind_metal = DM8transInd(iy, ix);
                    ind_diel  = DM9transInd(iy, ix);
                    FPMpoked(iy, ix) = mp.complexTransCompact(ind_diel, ind_metal, modvar.sbpIndex);
                end
            end            

            dEF3box = ((transOuterFPM-FPMpoked) - (transOuterFPM-fpm(eta_box_ind, xi_box_ind))) .* ...
                EF3inc(eta_box_ind, xi_box_ind); % Delta field (in a small region) at the FPM

            %--Matrices for the MFT from the FPM stamp to the Lyot stop
            rect_mat_pre = (exp(-2*pi*1j*(mp.P4.compact.ys*eta_box)/(lambda*mp.fl)))...
                *sqrt(mp.P4.compact.dx*mp.P4.compact.dx)*scaleFac*sqrt(mp.F3.compact.dxi*mp.F3.compact.deta)/(lambda*mp.fl);
            rect_mat_post  = (exp(-2*pi*1j*(xi_box*mp.P4.compact.xs)/(lambda*mp.fl)));

            %--MFT from FPM to Lyot stop (Nominal term transOuterFPM*EP4noFPM subtracts out to 0 since it ignores the FPM change).
            EP4 = 0 - rect_mat_pre*dEF3box*rect_mat_post; % MFT from FPM (F3) to Lyot stop plane (P4)
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering); %--Get the correct orientation 
            EP4 = mp.P4.compact.croppedMask .* EP4; %--Apply Lyot stop

            %--MFT to final focal plane
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.NrelayFend, mp.centering); %--Rotate the final image 180 degrees if necessary
            EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, mp.Fend.dxi, mp.Fend.Nxi, mp.Fend.deta, mp.Fend.Neta, mp.centering);
            if mp.useGPU; EFend = gather(EFend); end

            Gmode(:, Gindex) = mp.dm9.act_sens*(1/stepFac)*mp.dm9.weight*EFend(mp.Fend.corr.maskBool) /...
                sqrt(mp.Fend.compact.I00(modvar.sbpIndex));
        end
        Gindex = Gindex + 1;       
    end
end

if(mp.useGPU)
    Gmode = gather(Gmode);
end

end %--END OF FUNCTION
