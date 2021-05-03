% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Run the full-knowledge optical model and return the final E-field.
% - Not used by the estimator and controller.
% - Only used to create simulated intensity images.
%
% "Fourier" layout:  
% - For design or basic modeling only. 
% - Instead use a full model with Fresnel propagations for more realistic simulations.
%
% INPUTS
% ------
% mp : structure of model parameters
% lambda : wavelength in meters
% Ein : 2D input E-field at entrance
% normFac : intensity normalization factor 
%
%
% OUTPUTS
% -------
% Eout : 2-D electric field at final plane of optical layout
% varargout{1}==Efiber : E-field at final plane when a single mode fiber
% is used

function [Eout, varargout] = model_full_Fourier(mp, lambda, Ein, normFac)

mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.full.NdmPad;

if mp.flagRotation
    NrelayFactor = 1;
else
    NrelayFactor = 0; % zero out the number of relays
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--Change model values if the full model has a different value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isfield(mp, 'full'))
    if(isfield(mp.full, 'dm1'))
        if(isfield(mp.full.dm1, 'xc'));  mp.dm1.xc = mp.full.dm1.xc;  end % x-center location of DM1 surface [actuator widths]
        if(isfield(mp.full.dm1, 'yc'));  mp.dm1.yc = mp.full.dm1.yc;  end % y-center location of DM1 surface [actuator widths]
        if(isfield(mp.full.dm1, 'V0'));  mp.dm1.V = mp.dm1.V + mp.full.dm1.V0;  end % Add some extra starting command to the voltages  [volts]
    end
    if(isfield(mp.full, 'dm2'))
        if(isfield(mp.full.dm2, 'xc'));  mp.dm2.xc = mp.full.dm2.xc;  end % x-center location of DM2 surface [actuator widths]
        if(isfield(mp.full.dm2, 'yc'));  mp.dm2.yc = mp.full.dm2.yc;  end % y-center location of DM2 surface [actuator widths]
        if(isfield(mp.full.dm2, 'V0'));  mp.dm2.V = mp.dm2.V + mp.full.dm2.V0;  end % Add some extra starting command to the voltages  [volts]
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(any(mp.dm_ind==1));  DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.dx, NdmPad); else; DM1surf=zeros(NdmPad); end
if(any(mp.dm_ind==2));  DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.dx, NdmPad); else; DM2surf=zeros(NdmPad); end

pupil = pad_crop(mp.P1.full.mask, NdmPad);
Ein = pad_crop(Ein, NdmPad);

if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(mp.dm_ind==1)); DM1surf = gpuArray(DM1surf); end
    if(any(mp.dm_ind==2)); DM2surf = gpuArray(DM2surf); end
end

if(mp.flagDM1stop); DM1stop = pad_crop(mp.dm1.full.mask, NdmPad); else; DM1stop = 1; end
if(mp.flagDM2stop); DM2stop = pad_crop(mp.dm2.full.mask, NdmPad); else; DM2stop = 1; end

if(mp.flagDMwfe)
    if(any(mp.dm_ind==1));  Edm1WFE = exp(2*pi*1j/lambda.*pad_crop(mp.dm1.wfe, NdmPad, 'extrapval', 0)); else; Edm1WFE = ones(NdmPad); end
    if(any(mp.dm_ind==2));  Edm2WFE = exp(2*pi*1j/lambda.*pad_crop(mp.dm2.wfe, NdmPad, 'extrapval', 0)); else; Edm2WFE = ones(NdmPad); end
else
    Edm1WFE = ones(NdmPad);
    Edm2WFE = ones(NdmPad);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation from P1 to P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Start E-field at pupil plane P1
EP1 = pupil .* Ein;

%--Re-image to pupil plane P2
EP2 = propcustom_relay(EP1, NrelayFactor*mp.Nrelay1to2, mp.centering);

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
Edm1 = propcustom_PTP(EP2, mp.P2.full.dx*NdmPad, lambda, mp.d_P2_dm1);
Edm1 = Edm1WFE .* DM1stop .* exp(mirrorFac*2*pi*1j*DM1surf/lambda) .* Edm1;

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1, mp.P2.full.dx*NdmPad, lambda, mp.d_dm1_dm2); 
Edm2 = Edm2WFE .* DM2stop .* exp(mirrorFac*2*pi*1j*DM2surf/lambda) .* Edm2;

%--Back-propagate to effective pupil at P2
EP2eff = propcustom_PTP(Edm2, mp.P2.full.dx*NdmPad, lambda, -1*(mp.d_dm1_dm2 + mp.d_P2_dm1));

%--Re-image to pupil P3
EP3 = propcustom_relay(EP2eff, NrelayFactor*mp.Nrelay2to3, mp.centering);

%--Apply the apodizer mask (if there is one)
if mp.flagApod
    EP3 = mp.P3.full.mask .* pad_crop(EP3, mp.P3.full.Narr); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation from P3 to P4 depends on coronagraph type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
switch upper(mp.coro)
    case{'VORTEX', 'VC'}
        if mp.flagApod == false
            EP3 = pad_crop(EP3, 2^nextpow2(mp.P1.full.Narr)); %--Crop down if there isn't an apodizer mask
        end
        % Get FPM charge 
        if(numel(mp.F3.VortexCharge)==1)
            % single value indicates fully achromatic mask
            charge = mp.F3.VortexCharge;
        else
            % Passing an array for mp.F3.VortexCharge with
            % corresponding wavelengths mp.F3.VortexCharge_lambdas
            % represents a chromatic vortex FPM
            charge = interp1(mp.F3.VortexCharge_lambdas, mp.F3.VortexCharge, lambda, 'linear', 'extrap');
        end
        EP4 = propcustom_mft_Pup2Vortex2Pup(EP3, charge, mp.P1.full.Nbeam/2, 0.3, 5, mp.useGPU, mp.F3.VortexSpotDiam*(mp.lambda0/lambda), mp.F3.VortexSpotOffsets*(mp.lambda0/lambda));  %--MFTs
        % Undo the rotation inherent to propcustom_mft_Pup2Vortex2Pup.m
        if ~mp.flagRotation; EP4 = propcustom_relay(EP4, -1, mp.centering); end
        
        EP4 = pad_crop(EP4, mp.P4.full.Narr);

    case{'SPLC', 'FLC'}
        %--MFT from SP to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.full.dx, mp.F3.full.dxi, mp.F3.full.Nxi, mp.F3.full.deta, mp.F3.full.Neta, mp.centering); %--E-field incident upon the FPM
        EF3 = mp.F3.full.mask.*EF3inc; %--Apply FPM
        %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
        EP4 = propcustom_mft_FtoP(EF3, mp.fl, lambda, mp.F3.full.dxi, mp.F3.full.deta, mp.P4.full.dx, mp.P4.full.Narr, mp.centering); %--E-field incident upon the Lyot stop
        EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);
            
    case{'LC', 'APLC'}
        %--MFT from apodizer plane to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.full.dx, mp.F3.full.dxi, mp.F3.full.Nxi, mp.F3.full.deta, mp.F3.full.Neta, mp.centering);
        % Apply (1-FPM) for Babinet's principle later
        EF3 = (1-mp.F3.full.mask) .* EF3inc;
        % Use Babinet's principle at the Lyot plane. This is the term without the FPM.
        EP4noFPM = propcustom_relay(EP3, NrelayFactor*mp.Nrelay3to4, mp.centering); %--Propagate forward another pupil plane 
        EP4subtrahend = propcustom_mft_FtoP(EF3, mp.fl, lambda, mp.F3.full.dxi, mp.F3.full.deta, mp.P4.full.dx, mp.P4.full.Narr, mp.centering); 
        EP4subtrahend = propcustom_relay(EP4subtrahend, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);
        EP4 = pad_crop(EP4noFPM, mp.P4.full.Narr) - EP4subtrahend;

    case{'HLC'}
        
        switch mp.layout
            case{'fourier'}
                %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
                t_Ti_base = 0;
                t_Ni_vec = 0;
                t_PMGI_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
                pol = 2;
                [tCoef, ~] = falco_thin_film_material_def(lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, lambda*mp.FPM.d0fac, pol);
                transOuterFPM = tCoef;
            case{'fpm_scale'}
                transOuterFPM = mp.FPM.mask(1, 1); %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
                scaleFac = lambda/mp.lambda0; % Focal plane sampling varies with wavelength
            otherwise
                error('Invalid combination of mp.layout and mp.coro')
        end
        
        %--MFT from apodizer plane to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.full.dx, scaleFac*mp.F3.full.dxi, ...
            mp.F3.full.Nxi, scaleFac*mp.F3.full.deta, mp.F3.full.Neta, mp.centering);
        % Apply (1-FPM) for Babinet's principle later
        EF3 = (transOuterFPM-mp.FPM.mask).*EF3inc; %- transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
        % Use Babinet's principle at the Lyot plane.
        EP4noFPM = propcustom_relay(EP3, NrelayFactor*mp.Nrelay3to4, mp.centering); %--Propagate forward another pupil plane 
        EP4noFPM = transOuterFPM * pad_crop(EP4noFPM, mp.P4.full.Narr); %--Apply the phase and amplitude change from the FPM's outer complex transmission.
        EP4subtrahend = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.full.dxi, ...
            scaleFac*mp.F3.full.deta, mp.P4.full.dx, mp.P4.full.Narr, mp.centering);     
        EP4subtrahend = propcustom_relay(EP4subtrahend, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);
        EP4 = EP4noFPM - EP4subtrahend;

    case{'EHLC'}
        %--MFT from apodizer plane to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.full.dx, mp.F3.full.dxi, ...
            mp.F3.full.Nxi, mp.F3.full.deta, mp.F3.full.Neta, mp.centering);
        EF3 = mp.FPM.mask.*EF3inc; %--Apply FPM
        %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
        EP4 = propcustom_mft_FtoP(EF3, mp.fl, lambda, mp.F3.full.dxi, mp.F3.full.deta, ...
            mp.P4.full.dx, mp.P4.full.Narr, mp.centering); 
        EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);

    otherwise
        error('model_full_Fourier.m: Modely type\t %s\t not recognized.\n', mp.coro);
end

%--Remove the FPM completely if normalization value is being found in the vortex case
if normFac == 0 
    switch upper(mp.coro)
        case{'VORTEX', 'VC'}
            EP4 = propcustom_relay(EP3, NrelayFactor*mp.Nrelay3to4, mp.centering);
            EP4 = pad_crop(EP4, mp.P4.full.Narr);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Back to common propagation any coronagraph type   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Apply the Lyot stop
EP4 = mp.P4.full.croppedMask .* EP4;

%--MFT from Lyot Stop to final focal plane (i.e., P4 to Fend)
EP4 = propcustom_relay(EP4, NrelayFactor*mp.NrelayFend, mp.centering); %--Rotate the final image if necessary
EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.full.dx, mp.Fend.dxi, mp.Fend.Nxi, ...
    mp.Fend.deta, mp.Fend.Neta, mp.centering);

%--Don't apply FPM if normalization value is being found
if(normFac==0)
    Eout = EFend;  %--Don't normalize if normalization value is being found
else
    Eout = EFend/sqrt(normFac); %--Apply normalization
end

if mp.useGPU; Eout = gather(Eout); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Fiber-only propagation   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mp.flagFiber
    if mp.flagLenslet
        Efiber = cell(mp.Fend.Nlens, 1);
        sbpIndex = find(mp.sbp_centers == lambda);
        
        for nlens = 1:mp.Fend.Nlens
            EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, mp.Fend.dxi, mp.Fend.Nxi, ...
                mp.Fend.deta, mp.Fend.Neta, mp.centering, 'xfc', mp.Fend.x_lenslet_phys(nlens), 'yfc', mp.Fend.y_lenslet_phys(nlens));
            Elenslet = EFend.*mp.Fend.lenslet.mask;
            EF5 = propcustom_mft_PtoF(Elenslet, mp.lensletFL, lambda, mp.Fend.dxi, mp.F5.dxi, ...
                mp.F5.Nxi, mp.F5.deta, mp.F5.Neta, mp.centering);
            Efiber{nlens} = mp.F5.fiberMode(:, :, sbpIndex).*sum(sum(mp.F5.fiberMode(:, :, sbpIndex).*conj(EF5)));
        end
        
        Efiber = permute(reshape(cell2mat(Efiber)', mp.F5.Nxi, mp.F5.Neta, mp.Fend.Nlens), [2, 1, 3]);
        varargout{1} = Efiber;
        
    else  % Fibers placed in the focal plane with no lenslets
        EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, mp.Fend.dxi, ...
            mp.Fend.Nxi, mp.Fend.deta, mp.Fend.Neta, mp.centering);

        sbpIndex = find(mp.sbp_centers == lambda);
        
        Efiber = zeros(mp.Fend.Nxi, mp.Fend.Neta);
        for i=1:mp.Fend.Nfiber
            Eonefiber = mp.Fend.fiberMode(:, :, sbpIndex, i) .* ...
                sum(sum(mp.Fend.fiberMode(:, :, sbpIndex, i).*conj(EFend)));
            Efiber = Efiber + Eonefiber;
        end
        
        varargout{1} = Efiber;

        figure(901);
        imagesc(log10(abs(Efiber).^2)); axis equal tight; colorbar;
    end
end

end %--END OF FUNCTION
