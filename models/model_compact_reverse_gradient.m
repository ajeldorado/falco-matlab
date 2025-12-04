% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% 
% Compact reverse (ajoint) compact model used for algorithmic differentiation.
%
% In the "compact model,"
% aberrations are condensed to the input pupil (and possibly also the 
% output pupil) after phase retrieval, COFFEE, or similar since aberrations
% cannot accurately enough be attributed to individual optics. That makes 
% this model primarily a Fourier model, with Fresnel propagation only 
% useful and needed between the deformable mirrors and key coronagraph 
% optics known to be out of the pupil or focal planes.
%
% INPUTS
% ------
% command_vec : vectorized new delta commands to the DMs
% mp : structure of model parameters
% lambda : wavelength in meters
% EestAll : 1-D vectorized estimated E-fields at detector from last WFSC iteration
% EFendPrev : 2D previously computed E-field from the compact model before applying delta DMs
% log10reg : The log10() of the regularization to use when computing the control command.
%
% OUTPUTS
% -------
% total_cost : AD-EFC cost function value including terms for intensity and DM usage.
% gradient : The output command vector, consisting of the vectorized and 
%            concatenated commands to both DMs.
%
% NOTES
% -----
%

function [total_cost, gradient] = model_compact_reverse_gradient(...
    command_vec, log10reg, mp, EestAll, EFend_list, Edm1post_list, Edm2pre_list, DM2surf_list)

    mirrorFac = 2.;  % Phase change is twice the DM surface height.
    NdmPad = mp.compact.NdmPad;

    % Complex transmission of points outside FPM
    switch mp.coro
        case 'HLC'
            transOuterFPM = mp.F3.compact.mask(1, 1);
        otherwise
            transOuterFPM = 1;
    end

    if mp.flagRotation
        NrelayFactor = 1;
    else
        NrelayFactor = 0;  % zero out the number of relays
    end

    surf_dm1_bar_total = 0;
    surf_dm2_bar_total = 0;

    % Store copies of initial cumulative DM commands
    mp.dm1.V0 = mp.dm1.V;
    mp.dm2.V0 = mp.dm2.V;

    % Calculate new cumulative DM commands and put into 2-D arrays
    dDM1Vvec = zeros(mp.dm1.NactTotal, 1);
    dDM2Vvec = zeros(mp.dm2.NactTotal, 1);
    dDM1Vvec(mp.dm1.act_ele) = command_vec(mp.ctrl.uLegend == 1);
    dDM2Vvec(mp.dm2.act_ele) = command_vec(mp.ctrl.uLegend == 2);
    dv_dm1 = reshape(dDM1Vvec, mp.dm1.Nact, mp.dm1.Nact);
    dv_dm2 = reshape(dDM2Vvec, mp.dm2.Nact, mp.dm2.Nact);

    total_cost = 0;  % initialize

    for imode = 1:mp.jac.Nmode
    
        modvar = ModelVariables;
        modvar.whichSource = 'star';
        modvar.sbpIndex = mp.jac.sbp_inds(imode);
        modvar.zernIndex = mp.jac.zern_inds(imode);
        modvar.starIndex = mp.jac.star_inds(imode);

        wvl = mp.sbp_centers(modvar.sbpIndex);
        kk = mirrorFac*2*pi/wvl;
        I00 = mp.Fend.compact.I00(modvar.sbpIndex);
        EestVec = EestAll(:, imode);
        Eest2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
        Eest2D(mp.Fend.corr.maskBool) = EestVec;

        % Get model-based E-field before deltas. Should be pre-computed for speed.
        EFendA = EFend_list(:, imode);
        Edm1post = Edm1post_list(:, :, imode);
        Edm2pre = Edm2pre_list(:, :, imode);
        DM2surf = DM2surf_list(:, :, imode);

        % Get model-based E-field With delta DM commands applied.
        mp.dm1.V = mp.dm1.V0 + dv_dm1;
        mp.dm2.V = mp.dm2.V0 + dv_dm2;
        EFendB = compact(mp, modvar);
        % Reset DM commands
        mp.dm1.V = mp.dm1.V0;
        mp.dm2.V = mp.dm2.V0;

        % Compute the delta E-field from the latest commands (model new - model old).
        dEend = EFendB - EFendA;

        EdhNew = Eest2D + dEend;
        DH = EdhNew(mp.Fend.corr.maskBool);
        int_in_dh = sum(abs(DH(:)).^2);
        total_cost = total_cost + mp.jac.weights(imode) * int_in_dh;

        % Gradient
        Fend_masked = mp.jac.weights(imode)*2/sqrt(I00)*EdhNew.*mp.Fend.corr.maskBool;

        EP4_grad = propcustom_mft_FtoP(Fend_masked, -mp.fl, wvl, mp.Fend.dxi, mp.Fend.deta, mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering);
        EP4_grad = propcustom_relay(EP4_grad, NrelayFactor*mp.NrelayFend, mp.centering);
        EP4_grad = EP4_grad .* conj(pad_crop(mp.P4.compact.croppedMask, mp.P4.compact.Narr));

        switch upper(mp.coro)
        
            case{'VORTEX', 'VC', 'AVC'}
    
                EP4_grad = pad_crop(EP4_grad, NdmPad);
                % There is 1 rotation inherent to propcustom_mft_PtoFtoP
                EP4_grad = propcustom_relay(EP4_grad, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);
    
                % Get phase scale factor for the FPM. 
                if numel(mp.F3.phaseScaleFac) == 1
                    % Single value indicates fully achromatic mask
                    phaseScaleFac = mp.F3.phaseScaleFac;
                else
                    % Passing an array for mp.F3.phaseScaleFac with corresponding
                    % wavelengthsin mp.F3.phaseScaleFacLambdas represents a
                    % chromatic phase FPM.
                    phaseScaleFac = interp1(mp.F3.phaseScaleFacLambdas, mp.F3.phaseScaleFac, lambda, 'linear', 'extrap');
                end
               
                if mp.F3.flagDimple
                    inVal = mp.F3.inVal;
                    outVal = mp.F3.outVal;
                    pixPerLamD = mp.F3.compact.res;
        
                    inputs.type = 'sawtooth';
                    inputs.N = ceil_even(pixPerLamD*mp.P1.compact.Nbeam);
                    inputs.charge = mp.F3.VortexCharge;
                    inputs.phaseScaleFac = phaseScaleFac;
                    inputs.clocking = mp.F3.clocking;
                    inputs.roddierradius = mp.F3.roddierradius;
                    inputs.roddierphase = mp.F3.roddierphase;
                    
                    inputs.res = mp.F3.compact.res;
                    FPMcoarse = falco_gen_azimuthal_phase_mask(inputs);
                    
                    inputs.type = mp.F3.phaseMaskType;
                    inputs.res = floor(pixPerLamD*mp.P1.compact.Nbeam/(2*mp.F3.outVal));
                    FPMfine = falco_gen_azimuthal_phase_mask(inputs); clear inputs;
                    
                    EP3_grad = propcustom_mft_PtoFtoP_multispot(EP4_grad, FPMcoarse, FPMfine, mp.P1.compact.Nbeam/2, inVal, outVal, mp.useGPU, 'for_reverse_gradient');
                else
    
                    inVal = mp.F3.inVal;
                    outVal = mp.F3.outVal;
                    spotDiam = mp.F3.VortexSpotDiam * (mp.lambda0/lambda);
                    spotOffsets = mp.F3.VortexSpotOffsets * (mp.lambda0/lambda);
                    pixPerLamD = mp.F3.compact.res;
                    
                    inputs.type = mp.F3.phaseMaskType;
                    inputs.N = ceil_even(pixPerLamD*mp.P1.compact.Nbeam);
                    inputs.charge = mp.F3.VortexCharge;
                    inputs.phaseScaleFac = phaseScaleFac;
                    inputs.clocking = mp.F3.clocking;
                    inputs.Nsteps = mp.F3.NstepStaircase;
                    fpm = falco_gen_azimuthal_phase_mask(inputs); clear inputs;
    
                    flagReverseGradient = true;
                    EP3_grad = propcustom_mft_PtoFtoP(EP4_grad, fpm, mp.P1.compact.Nbeam/2, inVal, outVal, mp.useGPU, spotDiam, spotOffsets, flagReverseGradient);
                end

            case {'FLC', 'SPLC'}
    
                EP4_grad = propcustom_relay(EP4_grad, NrelayFactor*mp.Nrelay3to4-1, mp.centering);
                EF3_grad = propcustom_mft_PtoF(EP4_grad, -mp.fl, wvl, mp.P2.compact.dx, mp.F3.compact.dxi, mp.F3.compact.Nxi, mp.F3.compact.deta, mp.F3.compact.Neta, mp.centering);  % E-field incident upon the FPM
                EF3_grad = conj(mp.F3.compact.mask) .* EF3_grad;
                EP3_grad = propcustom_mft_FtoP(EF3_grad, -mp.fl, wvl, mp.F3.compact.dxi, mp.F3.compact.deta, mp.P4.compact.dx, NdmPad, mp.centering);  % Subtrahend term for Babinet's principle 
    
            case{'LC', 'APLC', 'HLC'}
                % MFT Method
                EP4noFPM_grad = EP4_grad;
                EP3noFPM_grad = pad_crop(propcustom_relay(EP4noFPM_grad, NrelayFactor*mp.Nrelay3to4, mp.centering), NdmPad);
                EP4subtr_grad = -EP4_grad;
                EP4subtr_grad = propcustom_relay(EP4subtr_grad, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);
                auxEP4subtr_grad = propcustom_mft_PtoF(EP4subtr_grad, -mp.fl, wvl, mp.P2.compact.dx, mp.F3.compact.dxi, mp.F3.compact.Nxi, mp.F3.compact.deta, mp.F3.compact.Neta, mp.centering);  % E-field incident upon the FPM
                EF3_grad = conj(transOuterFPM - mp.F3.compact.mask) .* auxEP4subtr_grad;
                EP3subtr_grad = propcustom_mft_FtoP(EF3_grad, -mp.fl, wvl, mp.F3.compact.dxi, mp.F3.compact.deta, mp.P4.compact.dx, NdmPad, mp.centering);  % Subtrahend term for Babinet's principle 
                EP3_grad = EP3noFPM_grad + EP3subtr_grad;
    
            otherwise
                error(sprintf('%s value of mp.coro not supported yet', mp.coro))

        end
        
        if mp.flagApod
            EP3_grad = conj(pad_crop(mp.P3.compact.mask, EP3_grad.shape)) .* EP3_grad;
        end

        EP2eff_grad = propcustom_relay(EP3_grad, NrelayFactor*mp.Nrelay2to3, mp.centering);

        % To DM2
        d_p2_to_dm2 = mp.d_P2_dm1 + mp.d_dm1_dm2;
        if d_p2_to_dm2 == 0
            Edm2_grad = EP2eff_grad;
        else
            Edm2_grad = falco.prop.ptp(EP2eff_grad, mp.P2.compact.dx*size(EP2eff_grad, 1), wvl, d_p2_to_dm2);
        end

        Edm2_grad = conj(exp(1j*kk*DM2surf)) * Edm2_grad;
        surf_dm2_bar = -kk*imag(Edm2_grad * conj(Edm2pre));

        % To DM1
        Edm1_grad = falco.prop.ptp(Edm2_grad, mp.P2.compact.dx*size(Edm2_grad, 1), wvl, -mp.d_dm1_dm2);
        surf_dm1_bar = -kk*imag(Edm1_grad * conj(Edm1post));

        surf_dm2_bar_total = surf_dm2_bar_total + mp.jac.weights(imode) * surf_dm2_bar;
        surf_dm1_bar_total = surf_dm1_bar_total + mp.jac.weights(imode) * surf_dm1_bar;
    
    end

    % Calculate DM penalty term component of cost function
    utu_coef = mp.ctrl.ad.utu_scale_fac * 10.0^(log10reg);
    total_cost = total_cost + utu_coef*sum(command_vec(:).^2);

    % TODO: Determine how to make the render_backprop method
    if mp.dm1.useDifferentiableModel && mp.dm2.useDifferentiableModel
        Vout1 = mp.dm1.differentiableModel.render_backprop(...
            surf_dm1_bar_total, wfe=False);
        Vout2 = mp.dm2.differentiableModel.render_backprop(...
            surf_dm2_bar_total, wfe=False);
    else
        error('mp.dm1.useDifferentiableModel and mp.dm2.useDifferentiableModel must be true for AD-EFC.')
    end
    Vout1 = Vout1 .* mp.dm1.VtoH;
    Vout2 = Vout2 .* mp.dm2.VtoH;
    gradient = [Vout1(mp.dm1.act_ele).'; Vout2(mp.dm2.act_ele).'];

end
