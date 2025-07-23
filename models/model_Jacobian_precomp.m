% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Precompute the starting point E-fields for all Jacobian modes.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% ------
% Gmode : Jacobian for the specified Zernike mode, DM, star, and subband.

function mp = model_Jacobian_precomp(mp)

    NdmPad = mp.compact.NdmPad;
    
    % % For vortex only
    % switch upper(mp.coro)
    %     case{'VORTEX', 'VC', 'AVC'}
    %         mp.jac.Nfft1_list = zeros(mp.jac.Nmode, 1);
    %         mp.jac.Nfft2_list = zeros(mp.jac.Nmode, 1);
    %         mp.jac.vortexDM1_list = zeros(Nfft1, Nfft1, mp.jac.Nmode);
    %         mp.jac.vortexDM2_list = zeros(Nfft1, Nfft1, mp.jac.Nmode);
    % end
    
    for imode = mp.jac.Nmode:-1:1
    
        modvar = ModelVariables;
        modvar.sbpIndex = mp.jac.sbp_inds(imode);
        modvar.zernIndex = mp.jac.zern_inds(imode);
        modvar.starIndex = mp.jac.star_inds(imode);
        lambda = mp.sbp_centers(modvar.sbpIndex); 
        surfIntoPhase = 2;
    
        % scaleFac = 1; % Default is that F3 focal plane sampling does not vary with wavelength
        
        if mp.flagRotation
            NrelayFactor = 1;
        else
            NrelayFactor = 0; % zero out the number of relays
        end

        % Define the FPM
        scaleFac = 1; % Default is that F3 focal plane sampling does not vary with wavelength
        switch upper(mp.coro)
            
            case{'LC', 'APLC', 'FLC', 'SPLC'}
                fpm = mp.F3.compact.mask;
                transOuterFPM = 1;
                
            case{'HLC'}
                switch mp.layout
                    case{'fourier'}
                        %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
                        if isfield(mp.compact, 'FPMcube')
                            fpm = squeeze(mp.compact.FPMcube(:, :, modvar.sbpIndex)); %--Complex transmission of the FPM.
                            transOuterFPM = fpm(1, 1);
                        else
                            t_Ti_base = 0;
                            t_Ni_vec = 0;
                            t_diel_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
                            pol = 2;
                            [transOuterFPM, ~] = falco_thin_film_material_def(mp.F3.substrate, mp.F3.metal, mp.F3.dielectric, lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_diel_vec, lambda*mp.FPM.d0fac, pol);
                            fpm = squeeze(mp.FPMcube(:, :, modvar.sbpIndex)); %--Complex transmission of the FPM. Calculated in model_Jacobian.m
                        end
                    case{'fpm_scale', 'proper', 'roman_phasec_proper', 'wfirst_phaseb_proper'}
                        fpm = squeeze(mp.compact.FPMcube(:, :, modvar.sbpIndex)); %--Complex transmission of the FPM.
                        transOuterFPM = fpm(1, 1); %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
                        scaleFac = lambda/mp.lambda0; % Focal plane sampling varies with wavelength
                    otherwise
                        transOuterFPM = 1;
                        fpm = 1;
                end

            otherwise
                transOuterFPM = 1;
                fpm = 1;
        end

        mp.jac.lc_fpm(:, :, imode) = fpm;
        mp.jac.transOuterFPM(imode) = transOuterFPM;
        mp.jac.scaleFac(imode) = scaleFac;
        
        switch upper(mp.coro)
            case{'VORTEX', 'VC', 'AVC'}
                %--Minimum FPM resolution for Jacobian calculations (in pixels per lambda/D)
                minPadFacVortex = 8; 
                
                %--Get phase scale factor for the FPM. 
                if numel(mp.F3.phaseScaleFac) == 1
                    % Single value indicates fully achromatic mask
                    phaseScaleFac = mp.F3.phaseScaleFac;
                else
                    % Passing an array for mp.F3.phaseScaleFac with corresponding
                    % wavelengthsin mp.F3.phaseScaleFacLambdas represents a
                    % chromatic phase FPM.
                    phaseScaleFac = interp1(mp.F3.phaseScaleFacLambdas, mp.F3.phaseScaleFac, lambda, 'linear', 'extrap');
                end
    
                %--Array size for planes P3, F3, and P4
                Nfft1 = 2^(ceil(log2(max([mp.dm1.compact.NdmPad, minPadFacVortex*mp.dm1.compact.Nbox])))); %--Don't crop--but do pad if necessary.
                mp.jac.Nfft1_list(imode) = Nfft1;
    
                %--Array size for planes P3, F3, and P4
                Nfft2 = 2^(ceil(log2(max([mp.dm2.compact.NdmPad, minPadFacVortex*mp.dm2.compact.Nbox])))); %--Don't crop--but do pad if necessary.
                mp.jac.Nfft2_list(imode) = Nfft2;
    
                % Variables related to the propagation method to/from the vortex mask
                NboxPad2AS = mp.dm2.compact.NboxAS;
                if mp.jac.mftToVortex
                    wBox = NboxPad2AS * mp.P2.compact.dx;
                    NactPerMat = wBox/mp.dm2.dm_spacing;
                    N = 6;
                    Nxi = ceil_even(NactPerMat*N*2);
                    pixres = N;
                    dxi = (mp.fl*lambda/wBox) / pixres;
                    deta = dxi;
                    if strcmpi(mp.centering, 'pixel')
                        xis = (-Nxi/2:(Nxi/2-1))*dxi;
                    else
                        xis = (-(Nxi-1)/2:(Nxi-1)/2)*dxi;
                    end
                    etas = xis.';

                    inputs.type = mp.F3.phaseMaskType;
                    inputs.N = Nxi;
                    inputs.charge = mp.F3.VortexCharge;
                    inputs.phaseScaleFac = phaseScaleFac;
                    inputs.clocking = mp.F3.clocking;
                    inputs.Nsteps = mp.F3.NstepStaircase;
                    inputs.res = pixres;
                    fpm = falco_gen_azimuthal_phase_mask(inputs); clear inputs;
            
                    mp.jac.vortexDM1_list(:, :, imode) = fpm;
                    mp.jac.vortexDM2_list(:, :, imode) = fpm;
                else        
                    % Generate central opaque spot
                    for ii = [1, 2]
    
                        if ii == 1
                            Nfft = Nfft1;
                        elseif ii == 2
                            Nfft = Nfft2;
                        end
    
                        if mp.F3.VortexSpotDiam > 0
                            inputs.pixresFPM = Nfft/mp.P1.compact.Nbeam; %--pixels per lambda/D
                            inputs.rhoInner = mp.F3.VortexSpotDiam/2*(mp.lambda0/lambda); % radius of inner FPM amplitude spot (in lambda_c/D)
                            inputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
                            inputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
                            inputs.centering = 'pixel';
                            spot = falco_gen_annular_FPM(inputs); clear inputs;
                            spot = pad_crop(spot, Nfft, 'extrapval', 1);
                        else
                            spot = 1;
                        end
                        
                        inputs.type = mp.F3.phaseMaskType;
                        inputs.N = Nfft;
                        inputs.charge = mp.F3.VortexCharge;
                        inputs.phaseScaleFac = phaseScaleFac;
                        inputs.clocking = mp.F3.clocking;
                        inputs.Nsteps = mp.F3.NstepStaircase;
                        inputs.res = Nfft/mp.P1.compact.Nbeam;
                        fpm = falco_gen_azimuthal_phase_mask(inputs); clear inputs;
                        
                        % Generate FPM with fftshift already applied
                        fftshiftVortex = fftshift(spot.*fpm);
    
                        if ii == 1
                            mp.jac.vortexDM1_list(:, :, imode) = fftshiftVortex;
                        elseif ii == 2
                            mp.jac.vortexDM2_list(:, :, imode) = fftshiftVortex;
                        end
    
                    end
                end
    
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
        
        % % Compute the change in E-field to apply at the exit pupil plane, EP4.
        % EP4mult = mp.P4.compact.E(:, :, modvar.sbpIndex);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Masks and DM surfaces
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        pupil = pad_crop(mp.P1.compact.mask, NdmPad);
        Ein = pad_crop(Ein, NdmPad);
        if(mp.useGPU)
            pupil = gpuArray(pupil);
            Ein = gpuArray(Ein);
        end
        
        % %--Re-image the apodizer from pupil P3 back to pupil P2.
        % if mp.flagApod
        %     apodReimaged = pad_crop(mp.P3.compact.mask, NdmPad);
        %     apodReimaged = propcustom_relay(apodReimaged, NrelayFactor*mp.Nrelay2to3, mp.centering);
        % else
        %     apodReimaged = ones(NdmPad); 
        % end
        
        if mp.flagDM1stop; DM1stop = pad_crop(mp.dm1.compact.mask, NdmPad); else; DM1stop = ones(NdmPad); end
        if mp.flagDM2stop; DM2stop = pad_crop(mp.dm2.compact.mask, NdmPad); else; DM2stop = ones(NdmPad); end
        mp.jac.DM2stop = pad_crop(DM2stop, mp.dm1.compact.NdmPad);
        
        if any(mp.dm_ind == 1); DM1surf = pad_crop(mp.dm1.compact.surfM, NdmPad);  else; DM1surf = 0; end 
        if any(mp.dm_ind == 2); DM2surf = pad_crop(mp.dm2.compact.surfM, NdmPad);  else; DM2surf = 0; end 
        if mp.useGPU
            if any(mp.dm_ind == 1); DM1surf = gpuArray(DM1surf); end
            if any(mp.dm_ind == 2); DM2surf = gpuArray(DM2surf); end
        end
        % mp.jac.DM1surf = DM1surf;
        mp.jac.DM2surf = pad_crop(mp.dm2.compact.surfM, mp.dm1.compact.NdmPad);
        
        % % Define the FPM
        % switch upper(mp.coro)
        % 
        %     case{'LC', 'APLC', 'FLC', 'SPLC'}
        %         fpm = mp.F3.compact.mask;
        %         transOuterFPM = 1;
        % 
        %     case{'HLC'}
        %         switch mp.layout
        %             case{'fourier'}
        %                 %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
        %                 if isfield(mp.compact, 'FPMcube')
        %                     fpm = squeeze(mp.compact.FPMcube(:, :, modvar.sbpIndex)); %--Complex transmission of the FPM.
        %                     transOuterFPM = fpm(1, 1);
        %                 else
        %                     t_Ti_base = 0;
        %                     t_Ni_vec = 0;
        %                     t_diel_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
        %                     pol = 2;
        %                     [transOuterFPM, ~] = falco_thin_film_material_def(mp.F3.substrate, mp.F3.metal, mp.F3.dielectric, lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_diel_vec, lambda*mp.FPM.d0fac, pol);
        %                     fpm = squeeze(mp.FPMcube(:, :, modvar.sbpIndex)); %--Complex transmission of the FPM. Calculated in model_Jacobian.m
        %                 end
        %             case{'fpm_scale', 'proper', 'roman_phasec_proper', 'wfirst_phaseb_proper'}
        %                 fpm = squeeze(mp.compact.FPMcube(:, :, modvar.sbpIndex)); %--Complex transmission of the FPM.
        %                 transOuterFPM = fpm(1, 1); %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
        %                 scaleFac = lambda/mp.lambda0; % Focal plane sampling varies with wavelength
        %             otherwise
        %                 error('Invalid combination of mp.layout and mp.coro')
        %         end
        % 
        %     otherwise
        %         error('Value of mp.coro not recognized.');
        % end
        
        % xisF3 = scaleFac * mp.F3.compact.xis;
        % etasF3 = scaleFac * mp.F3.compact.etas;
        % dxiF3 = scaleFac * mp.F3.compact.dxi;
        % detaF3 = scaleFac * mp.F3.compact.deta;
        
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
            % Store outputs
            mp.jac.Eocculted_list(:, :, imode) = Eocculted;
            mp.jac.Epeak_list(imode) = Epeak;

            if mp.P4.compact.Nbeam ~= mp.P1.compact.Nbeam
                if ~isfield(mp.P4.compact, 'maskAtP1res')
                    error('For peak Jacobian calculation, there must be a Lyot stop named mp.P4.compact.maskAtP1res that is sampled at the same resolution as the input pupil.')
                else
                    lyotStopReimaged = propcustom_relay(pad_crop(mp.P4.compact.maskAtP1res, NdmPad), NrelayFactor*(mp.Nrelay2to3+mp.Nrelay3to4));
                end
                
            else
                lyotStopReimaged = propcustom_relay(pad_crop(mp.P4.compact.mask, NdmPad), NrelayFactor*(mp.Nrelay2to3+mp.Nrelay3to4));
            end
            mp.jac.lyotStopReimaged = lyotStopReimaged;

        end
        
        %--Define pupil P1 and Propagate to pupil P2
        EP1 = pupil .* Ein; %--E-field at pupil plane P1
        EP2 = propcustom_relay(EP1, NrelayFactor*mp.Nrelay1to2, mp.centering); 
        
        %--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
        Edm1 = propcustom_PTP(EP2, mp.P2.compact.dx*NdmPad, lambda, mp.d_P2_dm1);
        Edm1 = Edm1WFE .* DM1stop .* exp(surfIntoPhase*2*pi*1j*DM1surf/lambda) .* Edm1;
    
        NboxPad2AS = mp.dm2.compact.NboxAS; 
        mp.dm2.compact.xy_box_lowerLeft_AS = mp.dm2.compact.xy_box_lowerLeft - (NboxPad2AS-mp.dm2.compact.Nbox)/2; %--Account for the padding of the influence function boxes
        
        % apodReimaged = pad_crop(apodReimaged, mp.dm2.compact.NdmPad);
        DM2stopPad = pad_crop(DM2stop, mp.dm2.compact.NdmPad);
        Edm2WFEpad = pad_crop(Edm2WFE, mp.dm2.compact.NdmPad);
        
        %--Propagate full field to DM2 before back-propagating in small boxes
        Edm2inc = propcustom_PTP(Edm1, mp.compact.NdmPad*mp.P2.compact.dx, lambda, mp.d_dm1_dm2); % E-field incident upon DM2
        Edm2 = pad_crop(Edm2inc, mp.dm2.compact.NdmPad);
        Edm2 = Edm2WFEpad.*DM2stopPad.*Edm2.*exp(surfIntoPhase*2*pi*1j/lambda*pad_crop(DM2surf, mp.dm2.compact.NdmPad)); % Initial E-field at DM2 including its own phase contribution
        
        % Store outputs
        mp.jac.Edm1_list(:, :, imode) = Edm1;
        mp.jac.Edm2_list(:, :, imode) = Edm2;
    end

end %--END OF FUNCTION
