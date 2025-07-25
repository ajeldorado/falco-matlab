% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% 
% Produce the electric field based on the knowledge
% available to the user in what is called the "compact model." Usually, 
% aberrations are condensed to the input pupil (and possibly also the 
% output pupil) after phase retrieval, COFFEE, or similar since aberrations
% cannot accurately enough be attributed to individual optics. That makes 
% this model primarily a Fourier model, with Fresnel propagation only 
% useful and needed between the deformable mirrors and key coronagraph 
% optics known to be out of the pupil or focal planes.
%
% INPUTS
% ------
% mp : structure of model parameters
% lambda : wavelength in meters
% Ein : 2D input E-field at entrance pupil
% normFac : intensity normalization factor 
% flagEval : boolean flag whether to evaluate at higher resolution for
%              throughput computation
% flagUseFPM : boolean flag whether to have the FPM in the beam or not
%
% OUTPUTS
% -------
% Eout : 2-D complex electric field in the final focal plane
%
% NOTES
% -----
% In wrapper above this that chooses layout, need to define mp.FPM.mask like this:
% mp.FPM.mask = falco_gen_HLC_FPM_complex_trans_mat(mp, modvar.sbpIndex, modvar.wpsbpIndex, 'compact');

function [Eout, Efiber, sDebug] = model_compact_general(mp, lambda, Ein, normFac, flagEval, flagUseFPM, varargin)

%--If there is an extra input, it is the exit pupil multiplier array.
EP4mult = 1; % default
if size(varargin, 2) == 1
    EP4mult = varargin{1};
end

if nargout >= 3
    debug = true;
else
    debug = false;
end

mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.compact.NdmPad;

switch mp.coro
    case 'HLC'
        transOuterFPM = mp.FPM.mask(1, 1);
    otherwise
        transOuterFPM = 1;
end

if mp.flagRotation
    NrelayFactor = 1;
else
    NrelayFactor = 0; % zero out the number of relays
end

if flagEval %--Higher resolution at final focal plane for computing stats such as throughput
    dxi = mp.Fend.eval.dxi;
    Nxi = mp.Fend.eval.Nxi;
    deta = mp.Fend.eval.deta;
    Neta = mp.Fend.eval.Neta; 
else %--Otherwise use the detector resolution
    dxi = mp.Fend.dxi;
    Nxi = mp.Fend.Nxi;
    deta = mp.Fend.deta;
    Neta = mp.Fend.Neta; 
end

% FPM scale factor with wavelength
switch mp.layout
    case 'fpm_scale'
        scaleFac = lambda/mp.lambda0; % Focal plane sampling varies with wavelength
    otherwise
        scaleFac = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Compute the DM surfaces for the current DM commands
if any(mp.dm_ind == 1); DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, NdmPad); else; DM1surf = 0; end
if any(mp.dm_ind == 2); DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, NdmPad); else; DM2surf = 0; end

pupil = pad_crop(mp.P1.compact.mask, NdmPad);
Ein = pad_crop(Ein, mp.compact.NdmPad);

if(mp.flagDM1stop); DM1stop = pad_crop(mp.dm1.compact.mask, NdmPad); else; DM1stop = ones(NdmPad); end
if(mp.flagDM2stop); DM2stop = pad_crop(mp.dm2.compact.mask, NdmPad); else; DM2stop = ones(NdmPad); end

if mp.useGPU
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if any(mp.dm_ind == 1); DM1surf = gpuArray(DM1surf); end
    if any(mp.dm_ind == 2); DM2surf = gpuArray(DM2surf); end
end

%--Apply WFE to DMs 1 and 2
if(mp.flagDMwfe)
    if any(mp.dm_ind == 1); Edm1WFE = exp(2*pi*1j/lambda.*pad_crop(mp.dm1.compact.wfe, NdmPad, 'extrapval', 0)); else; Edm1WFE = ones(NdmPad); end
    if any(mp.dm_ind == 2); Edm2WFE = exp(2*pi*1j/lambda.*pad_crop(mp.dm2.compact.wfe, NdmPad, 'extrapval', 0)); else; Edm2WFE = ones(NdmPad); end
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
if debug, sDebug.EP2_before_dms = EP2; end

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
Edm1 = propcustom_PTP(EP2, mp.P2.compact.dx*NdmPad, lambda, mp.d_P2_dm1);
Edm1 = Edm1WFE .* DM1stop .* exp(mirrorFac*2*pi*1j*DM1surf/lambda) .* Edm1;

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1, mp.P2.compact.dx*NdmPad, lambda, mp.d_dm1_dm2); 
Edm2 = Edm2WFE .* DM2stop .* exp(mirrorFac*2*pi*1j*DM2surf/lambda) .* Edm2;

%--Back-propagate to effective pupil at P2
EP2eff = propcustom_PTP(Edm2, mp.P2.compact.dx*NdmPad, lambda, -1*(mp.d_dm1_dm2 + mp.d_P2_dm1));
if debug, sDebug.EP2_after_dms = EP2eff; end

%--Re-image to pupil P3
EP3 = propcustom_relay(EP2eff, NrelayFactor*mp.Nrelay2to3, mp.centering);

%--Apply apodizer mask.
if mp.flagApod
    EP3 = mp.P3.compact.mask .* pad_crop(EP3, mp.P3.compact.Narr); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation from P3 to P4 depends on coronagraph type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Remove FPM from beam to get normalization constant for vortex coronagraphs
switch upper(mp.coro)
    case{'VORTEX', 'VC', 'AVC'}
        if normFac == 0
            flagUseFPM = false;
        end
end

if flagUseFPM
    switch upper(mp.coro)

        case{'VORTEX', 'VC'}

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
                
                EP4 = propcustom_mft_PtoFtoP_multispot(EP3, FPMcoarse, FPMfine, mp.P1.compact.Nbeam/2, inVal, outVal, mp.useGPU);
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

                EP4 = propcustom_mft_PtoFtoP(EP3, fpm, mp.P1.compact.Nbeam/2, inVal, outVal, mp.useGPU, spotDiam, spotOffsets);
            end

            % One 180-degree rotation is inherent to propcustom_mft_PtoFtoP
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);

            % Resize beam if Lyot plane has different resolution
            if mp.P4.compact.Nbeam ~= mp.P1.compact.Nbeam
                N1 = length(EP4);
                N4 = ceil_even(1.1*mp.P4.compact.Narr);
                if strcmpi(mp.centering, 'pixel')
                    x1 = (-N1/2:(N1/2-1))/mp.P1.compact.Nbeam;
                    x4 = (-N4/2:(N4/2-1))/mp.P4.compact.Nbeam;
                else
                    x1 = (-(N1-1)/2:(N1-1)/2)/mp.P1.compact.Nbeam;
                    x4 = (-(N4-1)/2:(N4-1)/2)/mp.P4.compact.Nbeam;
                end
                [X1, Y1] = meshgrid(x1);
                [X4, Y4] = meshgrid(x4);
                EP4 = interp2(X1, Y1, EP4, X4, Y4);
                % Preserve summed intensity in the pupil:
                EP4 = (mp.P1.compact.Nbeam/mp.P4.compact.Nbeam) * EP4;
            end

            EP4 = pad_crop(EP4, mp.P4.compact.Narr);

        case{'SPLC', 'FLC'}
            %--MFT from SP to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.compact.dx, scaleFac*mp.F3.compact.dxi,...
                mp.F3.compact.Nxi, scaleFac*mp.F3.compact.deta, mp.F3.compact.Neta, mp.centering);
            EF3 = mp.F3.compact.mask .* EF3inc; % Apply FPM
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4 = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.compact.dxi, scaleFac*mp.F3.compact.deta, mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering);
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);

        case{'LC', 'APLC'}
            %--MFT from SP to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.compact.dx, scaleFac*mp.F3.compact.dxi, ...
                mp.F3.compact.Nxi, scaleFac*mp.F3.compact.deta, mp.F3.compact.Neta, mp.centering);
            %--Apply (1-FPM) for Babinet's principle later
            EF3 = (1 - mp.F3.compact.mask) .* EF3inc;

            %--Use Babinet's principle at the Lyot plane.
            EP4noFPM = propcustom_relay(EP3, NrelayFactor*mp.Nrelay3to4, mp.centering); %--Propagate forward another pupil plane 
            EP4noFPM = pad_crop(EP4noFPM, mp.P4.compact.Narr); %--Crop down to the size of the Lyot stop opening
            EP4subtrahend = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.compact.dxi, ...
                scaleFac*mp.F3.compact.deta, mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering);
            EP4subtrahend = propcustom_relay(EP4subtrahend, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);
            EP4 = EP4noFPM - EP4subtrahend;

        case{'HLC'}
            
            switch mp.layout
                case{'fourier'}
                    
                    scaleFac = 1; % Focal plane sampling does not vary with wavelength
                    
                case{'fpm_scale', 'proper', 'roman_phasec_proper', 'wfirst_phaseb_proper'}
                    
                    scaleFac = lambda/mp.lambda0; % Focal plane sampling varies with wavelength
                    
                otherwise
                    
                    error('Invalid combination of mp.layout and mp.coro')
                    
            end

            %--Propagate to focal plane F3
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.compact.dx, scaleFac*mp.F3.compact.dxi,...
                mp.F3.compact.Nxi, scaleFac*mp.F3.compact.deta, mp.F3.compact.Neta, mp.centering);
            %--Apply (1-FPM) for Babinet's principle later
            EF3 = (transOuterFPM - mp.FPM.mask).*EF3inc; 
            %--Use Babinet's principle at the Lyot plane.
            EP4noFPM = propcustom_relay(EP3, NrelayFactor*mp.Nrelay3to4, mp.centering);
            EP4noFPM = transOuterFPM * pad_crop(EP4noFPM, mp.P4.compact.Narr); %--Apply the phase and amplitude change from the FPM's outer complex transmission.
            EP4subtrahend = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.compact.dxi, ...
                scaleFac*mp.F3.compact.deta, mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering);  
            EP4subtrahend = propcustom_relay(EP4subtrahend, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);
            EP4 = (EP4noFPM-EP4subtrahend); 

        case{'EHLC'}
            %--MFT from apodizer plane to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl, lambda, mp.P2.compact.dx, scaleFac*mp.F3.compact.dxi, ...
                mp.F3.compact.Nxi, scaleFac*mp.F3.compact.deta, mp.F3.compact.Neta, mp.centering);
            EF3 = mp.FPM.mask .* EF3inc; %--Apply FPM
            EP4 = propcustom_mft_FtoP(EF3, mp.fl, lambda, scaleFac*mp.F3.compact.dxi, scaleFac*mp.F3.compact.deta, ...
                mp.P4.compact.dx, mp.P4.compact.Narr, mp.centering);
            EP4 = propcustom_relay(EP4, NrelayFactor*mp.Nrelay3to4 - 1, mp.centering);

    end

else % No FPM in beam path, so relay directly from P3 to P4.
    
    EP4 = propcustom_relay(EP3, NrelayFactor*mp.Nrelay3to4, mp.centering);
    EP4 = transOuterFPM * EP4;

    % Interpolate beam if Lyot plane has different resolution
    if mp.P4.compact.Nbeam ~= mp.P1.compact.Nbeam
        %Make sure array is oversized before downsampling
        padFac = 1.2;
        EP4 = pad_crop(EP4, ceil_even(padFac*mp.P1.compact.Nbeam));
        N1 = length(EP4);
        N4 = ceil_even(padFac*mp.P4.compact.Narr);
        if strcmpi(mp.centering, 'pixel')
            x1 = (-N1/2:(N1/2-1)) / mp.P1.compact.Nbeam;
            x4 = (-N4/2:(N4/2-1)) / mp.P4.compact.Nbeam;
        elseif strcmpi(mp.centering, 'interpixel')
            x1 = (-(N1-1)/2:(N1-1)/2) / mp.P1.compact.Nbeam;
            x4 = (-(N4-1)/2:(N4-1)/2) / mp.P4.compact.Nbeam;
        end
        [X1, Y1] = meshgrid(x1);
        [X4, Y4] = meshgrid(x4);
        EP4 = interp2(X1, Y1, EP4, X4, Y4);
        % Preserve summed intensity in the pupil:
        EP4 = (mp.P1.compact.Nbeam/mp.P4.compact.Nbeam) * EP4;
    end
    
    EP4 = pad_crop(EP4, mp.P4.compact.Narr);
    
end

if debug, sDebug.EP4_before_mask = EP4; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Back to common propagation any coronagraph type   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Apply the Lyot stop
EP4 = mp.P4.compact.croppedMask .* EP4;

if debug, sDebug.EP4_after_mask = EP4; end

% Apply rotation, downstream pointing, and dowstream aberrations.
EP4 = EP4mult .* EP4;
EP4 = propcustom_relay(EP4, NrelayFactor*mp.NrelayFend, mp.centering); %--Rotate the final image if necessary


%--MFT to camera
EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, dxi, Nxi, deta, Neta, mp.centering);

%--Don't apply FPM if normalization value is being found
if normFac == 0
    Eout = EFend; %--Don't normalize if normalization value is being found
else
    Eout = EFend/sqrt(normFac); %--Apply normalization
end

if mp.useGPU
    Eout = gather(Eout);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Fiber-only propagation   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Efiber = 0;
if mp.flagFiber && ~flagEval
    if mp.flagLenslet % **not up to date with lenslets**
        Efiber = cell(mp.Fend.Nlens, 1);
        sbpIndex = find(mp.sbp_centers == lambda);
        
        for nlens = 1:mp.Fend.Nlens
            EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, mp.Fend.dxi, mp.Fend.Nxi, mp.Fend.deta, mp.Fend.Neta, mp.centering, 'xfc', mp.Fend.x_lenslet_phys(nlens), 'yfc', mp.Fend.y_lenslet_phys(nlens));
            Elenslet = EFend.*mp.Fend.lenslet.mask;
            EF5 = propcustom_mft_PtoF(Elenslet, mp.lensletFL, lambda, mp.Fend.dxi, mp.F5.dxi, mp.F5.Nxi, mp.F5.deta, mp.F5.Neta, mp.centering);
            Efiber{nlens} = mp.F5.fiberMode(:, :, sbpIndex).*sum(sum(mp.F5.fiberMode(:, :, sbpIndex).*conj(EF5)));
        end
        
        Efiber = permute(reshape(cell2mat(Efiber)', mp.F5.Nxi, mp.F5.Neta, mp.Fend.Nlens), [2, 1, 3]);
        varargout{1} = Efiber;
        
    else  %Fibers placed in the focal plane with no lenslets
%         EFend = propcustom_mft_PtoF(EP4, mp.fl, lambda, mp.P4.compact.dx, mp.Fend.dxi, mp.Fend.Nxi, mp.Fend.deta, mp.Fend.Neta, mp.centering);

        sbpIndex = find(mp.sbp_centers == lambda);
        
        Efiber = zeros(mp.Fend.Nfiber, 1);
        for ii=1:mp.Fend.Nfiber
%             if flagEval
%                 normFactor_fiber = sqrt(mp.Fend.eval.I00_fiber(ii,sbpIndex));
%             else
            normFactor_fiber = sqrt(mp.Fend.compact.I00_fiber(ii,sbpIndex));
%             end
            Eonefiber = sum(sum(mp.Fend.fiberMode(:, :, sbpIndex, ii).*EFend)) / normFactor_fiber;
            Efiber(ii) = Eonefiber;
        end
        varargout{1} = Efiber;
    end
end

end
