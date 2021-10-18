% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Wrapper for the simplified optical models used for the fast Jacobian calculation.
% The first-order derivative of the DM pokes are propagated through the system.
% Does not include unknown aberrations/errors that are in the full model.
%
% For computing the control Jacobian of the coronagraph with the FPM
% removed (but Lyot stop still in place).
%
% INPUTS
% ------
% mp : structure of model parameters
% iMode : index of the pair of sub-bandpass index and Zernike mode index
% whichDM : which DM number
%
% OUTPUTS
% -------
% Gmode : Jacobian for the specified Zernike mode, DM, star, and subband.

function Gmode = model_Jacobian_no_FPM(mp, iMode, whichDM)

modvar.sbpIndex = mp.jac.sbp_inds(iMode);
modvar.zernIndex = mp.jac.zern_inds(iMode);
modvar.starIndex = mp.jac.star_inds(iMode);
lambda = mp.sbp_centers(modvar.sbpIndex); 
NdmPad = mp.compact.NdmPad;
surfIntoPhase = 2;

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

EttUndoAtP2 = 1./propcustom_relay(Ett, NrelayFactor*mp.Nrelay1to2);

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

if mp.P4.compact.Nbeam ~= mp.P1.compact.Nbeam
    if ~isfield(mp.P4.compact, 'maskAtP1res')
        error(['For peak Jacobian calculation, there must be a Lyot stop named mp.P4.compact.maskAtP1res'
            'that is sampled at the same resolution as the input pupil.'])
    else
        lyotStopReimaged = propcustom_relay(pad_crop(mp.P4.compact.maskAtP1res, NdmPad), NrelayFactor*(mp.Nrelay2to3+mp.Nrelay3to4));
    end
    
else
    lyotStopReimaged = propcustom_relay(pad_crop(mp.P4.compact.mask, NdmPad), NrelayFactor*(mp.Nrelay2to3+mp.Nrelay3to4));
end

pupil = pad_crop(mp.P1.compact.mask, NdmPad);
Ein = pad_crop(Ein, NdmPad);

if mp.flagDM1stop; DM1stop = pad_crop(mp.dm1.compact.mask, NdmPad); else; DM1stop = ones(NdmPad); end
if mp.flagDM2stop; DM2stop = pad_crop(mp.dm2.compact.mask, NdmPad); else; DM2stop = ones(NdmPad); end

if any(mp.dm_ind == 1); DM1surf = pad_crop(mp.dm1.compact.surfM, NdmPad);  else; DM1surf = 0; end 
if any(mp.dm_ind == 2); DM2surf = pad_crop(mp.dm2.compact.surfM, NdmPad);  else; DM2surf = 0; end 

if mp.useGPU
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if any(mp.dm_ind == 1); DM1surf = gpuArray(DM1surf); end
    if any(mp.dm_ind == 2); DM2surf = gpuArray(DM2surf); end
end

%--Re-image the apodizer from pupil P3 back to pupil P2.
if mp.flagApod
    apodReimaged = pad_crop(mp.P3.compact.mask, NdmPad);
    apodReimaged = propcustom_relay(apodReimaged, -NrelayFactor*mp.Nrelay2to3, mp.centering);
else
    apodReimaged = ones(NdmPad); 
end

% Define the complex transmission of the outer part of the FPM
switch upper(mp.coro)
        
    case{'HLC'}        
        %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
        t_Ti_base = 0;
        t_Ni_vec = 0;
        t_PMGI_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
        pol = 2;
        [transOuterFPM, ~] = falco_thin_film_material_def(lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, lambda*mp.FPM.d0fac, pol);
        
    otherwise
        transOuterFPM = 1;
end

scaleFac = 1; % Default is that F3 focal plane sampling does not vary with wavelength
switch upper(mp.coro)
    case{'HLC'}
        switch mp.layout
            case{'fpm_scale', 'proper', 'roman_phasec_proper', 'wfirst_phaseb_proper'}
                scaleFac = lambda/mp.lambda0; % Focal plane sampling varies with wavelength
        end
end


%--For including DM surface errors (quilting, scalloping, etc.)
if mp.flagDMwfe
    if any(mp.dm_ind == 1); Edm1WFE = exp(2*pi*1j/lambda.*pad_crop(mp.dm1.compact.wfe, NdmPad, 'extrapval', 0)); else; Edm1WFE = ones(NdmPad); end
    if any(mp.dm_ind == 2); Edm2WFE = exp(2*pi*1j/lambda.*pad_crop(mp.dm2.compact.wfe, NdmPad, 'extrapval', 0)); else; Edm2WFE = ones(NdmPad); end
else
    Edm1WFE = ones(NdmPad);
    Edm2WFE = ones(NdmPad);
end

% % For interpolating beam if Lyot plane has different resolution
% if mp.P4.compact.Nbeam ~= mp.P1.compact.Nbeam
%     N1 = mp.P1.compact.Nbeam; %length(EP4);
%     % mag = mp.P4.compact.Nbeam / mp.P1.compact.Nbeam;
%     N4 = mp.P4.compact.Narr; %ceil_even(mag*N1);
%     if strcmpi(mp.centering, 'pixel')
%         x1 = (-N1/2:(N1/2-1)) / mp.P1.compact.Nbeam;
%         x4 = (-N4/2:(N4/2-1)) / mp.P4.compact.Nbeam;
%     elseif strcmpi(mp.centering, 'interpixel')
%         x1 = (-(N1-1)/2:(N1-1)/2) / mp.P1.compact.Nbeam;
%         x4 = (-(N4-1)/2:(N4-1)/2) / mp.P4.compact.Nbeam;
%     end
%     [X1, Y1] = meshgrid(x1);
%     [X4, Y4] = meshgrid(x4);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Get the unocculted peak E-field and coronagraphic E-field
% if mp.jac.minimizeNI
%     modvar.whichSource = 'star';
%     Eunocculted = model_compact(mp, modvar, 'nofpm');
%     [~, indPeak] = max(abs(Eunocculted(:)));
% end

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil .* Ein;
EP2 = propcustom_relay(EP1, NrelayFactor*mp.Nrelay1to2, mp.centering);

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
Edm1 = propcustom_PTP(EP2, mp.P2.compact.dx*NdmPad, lambda, mp.d_P2_dm1);
Edm1 = Edm1WFE .* DM1stop .* exp(surfIntoPhase*2*pi*1j*DM1surf/lambda) .* Edm1;

%--DM1---------------------------------------------------------
if whichDM == 1 
    Gmode = zeros(1, mp.dm1.Nele);
    
    %--Two array sizes (at same resolution) of influence functions for MFT and angular spectrum
    NboxPad1AS = mp.dm1.compact.NboxAS; %NboxPad1;%2.^ceil(log2(NboxPad1)); %--Power of 2 array size for FFT-AS propagations from DM1->DM2->DM1
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
    
    Edm1pad = pad_crop(Edm1, mp.dm1.compact.NdmPad); %--Pad or crop for expected sub-array indexing
    Edm2WFEpad = pad_crop(Edm2WFE, mp.dm1.compact.NdmPad); %--Pad or crop for expected sub-array indexing
    apodReimaged = pad_crop(apodReimaged, mp.dm1.compact.NdmPad);
    lyotStopReimaged = pad_crop(lyotStopReimaged, mp.dm1.compact.NdmPad);
    EttUndoAtP2 = pad_crop(EttUndoAtP2, mp.dm1.compact.NdmPad);

    %--Propagate each actuator from DM1 through the optical system
    Gindex = 1; % initialize index counter
    for iact = mp.dm1.act_ele(:).'  %--MUST BE A COLUMN VECTOR
        
        if any(any(mp.dm1.compact.inf_datacube(:, :, iact)))
            %--influence function indices of the padded influence function in the full padded pupil
            x_box_AS_ind = mp.dm1.compact.xy_box_lowerLeft_AS(1, iact):...
                mp.dm1.compact.xy_box_lowerLeft_AS(1, iact)+NboxPad1AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = mp.dm1.compact.xy_box_lowerLeft_AS(2, iact):...
                mp.dm1.compact.xy_box_lowerLeft_AS(2, iact)+NboxPad1AS-1; % y-indices in pupil arrays for the box
            
            %--Propagate from DM1 to DM2, and then back to P2
            dEbox = (surfIntoPhase*2*pi*1j/lambda)*pad_crop(mp.dm1.VtoH(iact)*mp.dm1.compact.inf_datacube(:, :, iact), NboxPad1AS); %--Pad influence function at DM1 for angular spectrum propagation.
            dEbox = propcustom_PTP_inf_func(dEbox.*Edm1pad(y_box_AS_ind, x_box_AS_ind), mp.P2.compact.dx*NboxPad1AS, lambda, mp.d_dm1_dm2, mp.dm1.dm_spacing, mp.propMethodPTP); % forward propagate to DM2 and apply DM2 E-field
            dEP2box = propcustom_PTP_inf_func(dEbox.*Edm2WFEpad(y_box_AS_ind, x_box_AS_ind).*DM2stop(y_box_AS_ind, x_box_AS_ind).*exp(surfIntoPhase*2*pi*1j/lambda*DM2surf(y_box_AS_ind, x_box_AS_ind)), mp.P2.compact.dx*NboxPad1AS, lambda, -1*(mp.d_dm1_dm2 + mp.d_P2_dm1), mp.dm1.dm_spacing, mp.propMethodPTP ); % back-propagate to DM1
            
            % Apply the reimaged apodizer at P2
            dEP2box = apodReimaged(y_box_AS_ind, x_box_AS_ind) .* dEP2box;
            
            % Put the star back on-axis (if it isn't already)
            dEP2box = EttUndoAtP2(y_box_AS_ind, x_box_AS_ind) .* dEP2box;
            
            % Apply the reimaged Lyot stop at P2
            dEP2box = lyotStopReimaged(y_box_AS_ind, x_box_AS_ind) .* dEP2box;
            
            dEFendPeak = sum(dEP2box(:)) * transOuterFPM *...
                sqrt(mp.P2.compact.dx*mp.P2.compact.dx) * ...
                scaleFac * sqrt(mp.Fend.dxi*mp.Fend.deta) / (lambda*mp.fl);
            
            Gmode(:, Gindex) = dEFendPeak / sqrt(mp.Fend.compact.I00(modvar.sbpIndex));
        end
        
        Gindex = Gindex + 1;
    end    

%--DM2---------------------------------------------------------
elseif whichDM == 2
    Gmode = zeros(1, mp.dm2.Nele);
    
    %--Two array sizes (at same resolution) of influence functions for MFT and angular spectrum
    NboxPad2AS = mp.dm2.compact.NboxAS; 
    mp.dm2.compact.xy_box_lowerLeft_AS = mp.dm2.compact.xy_box_lowerLeft - (NboxPad2AS-mp.dm2.compact.Nbox)/2; %--Account for the padding of the influence function boxes
    
    DM2stopPad = pad_crop(DM2stop, mp.dm2.compact.NdmPad);
    Edm2WFEpad = pad_crop(Edm2WFE, mp.dm2.compact.NdmPad);
    apodReimaged = pad_crop(apodReimaged, mp.dm2.compact.NdmPad);
    lyotStopReimaged = pad_crop(lyotStopReimaged, mp.dm2.compact.NdmPad);
    EttUndoAtP2 = pad_crop(EttUndoAtP2, mp.dm2.compact.NdmPad);
    
    %--Propagate full field to DM2 before back-propagating in small boxes
    Edm2inc = pad_crop(propcustom_PTP(Edm1, mp.compact.NdmPad*mp.P2.compact.dx, lambda, mp.d_dm1_dm2), mp.dm2.compact.NdmPad); % E-field incident upon DM2
    Edm2inc = pad_crop(Edm2inc, mp.dm2.compact.NdmPad);
    Edm2 = Edm2WFEpad.*DM2stopPad.*Edm2inc.*exp(surfIntoPhase*2*pi*1j/lambda*pad_crop(DM2surf, mp.dm2.compact.NdmPad)); % Initial E-field at DM2 including its own phase contribution
    
    %--Propagate each actuator from DM2 through the rest of the optical system
    Gindex = 1; % Initialize index counter
    for iact = mp.dm2.act_ele(:).'  %--Only compute for acutators specified %--MUST BE A COLUMN VECTOR
        
        if any(any(mp.dm2.compact.inf_datacube(:, :, iact)))
            %--influence function indices of the padded influence function in the full padded pupil
            x_box_AS_ind = mp.dm2.compact.xy_box_lowerLeft_AS(1, iact):mp.dm2.compact.xy_box_lowerLeft_AS(1, iact)+NboxPad2AS-1; % x-indices in pupil arrays for the box
            y_box_AS_ind = mp.dm2.compact.xy_box_lowerLeft_AS(2, iact):mp.dm2.compact.xy_box_lowerLeft_AS(2, iact)+NboxPad2AS-1; % y-indices in pupil arrays for the box

            dEbox = mp.dm2.VtoH(iact)*(surfIntoPhase*2*pi*1j/lambda)*pad_crop(mp.dm2.compact.inf_datacube(:, :, iact), NboxPad2AS); %--the padded influence function at DM2
            dEP2box = propcustom_PTP_inf_func(dEbox.*Edm2(y_box_AS_ind, x_box_AS_ind), mp.P2.compact.dx*NboxPad2AS, lambda, -1*(mp.d_dm1_dm2 + mp.d_P2_dm1), mp.dm2.dm_spacing, mp.propMethodPTP); % back-propagate to pupil P2

            % Apply the reimaged apodizer at P2
            dEP2box = apodReimaged(y_box_AS_ind, x_box_AS_ind) .* dEP2box;
            
            % Put the star back on-axis (if it isn't already)
            dEP2box = EttUndoAtP2(y_box_AS_ind, x_box_AS_ind) .* dEP2box;
            
            % Apply the reimaged Lyot stop at P2
            dEP2box = lyotStopReimaged(y_box_AS_ind, x_box_AS_ind) .* dEP2box;
            
            dEFendPeak = sum(dEP2box(:)) * transOuterFPM *...
                sqrt(mp.P2.compact.dx*mp.P2.compact.dx) * ...
                scaleFac * sqrt(mp.Fend.dxi*mp.Fend.deta) / (lambda*mp.fl);
            
            Gmode(:, Gindex) = dEFendPeak / sqrt(mp.Fend.compact.I00(modvar.sbpIndex));
        end
        
        Gindex = Gindex + 1;
    end

elseif whichDM == 8
    Gmode = zeros(1, mp.dm8.Nele);
    
elseif whichDM == 9
    Gmode = zeros(1, mp.dm8.Nele);
    
end

if mp.useGPU
    Gmode = gather(Gmode);
end

end %--END OF FUNCTION
