function Eout = model_ZWFS(mp, modvar, varargin)
% function Eout = model_ZWFS(mp, modvar, varargin)
% Function to model a Zernike WFS. 
% Options: 
%   - Default: propagates through primary, apodizer, DMs, Zernike WFS focal
%   plane mask, and returns the output pupil field.
%   - 'refwave' computes the reference wave, which only propagates the part
%       of the beam that passes through the phase dimple.
%   - 'to_input' only propagates to the input pupil.

% set default options 
refwave = false;
to_input = false;

% Check if user wants to compute the reference wave, which is the FT of the
% part of the beam that passes through the phase dimple. 
if(~isempty(varargin))
    flag = varargin(1); 
    if(strcmp(flag,'refwave'))
        refwave = true; 
    end
    if(strcmp(flag,'to_input'))
        to_input = true;
        refwave = false; 
    end
end

%--Set the wavelength
if(isfield(modvar,'lambda'))
    lambda = modvar.lambda;
elseif(isfield(modvar,'ebpIndex'))
    lambda = mp.wfs.lambdas(modvar.ebpIndex);
elseif(isfield(modvar,'sbpIndex'))
    lambda = mp.wfs.lambdasMat(modvar.sbpIndex,modvar.wpsbpIndex);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input E-fields
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For the out-of-band case, scale the E-field defined for the coronagraph 
% full model assuming the phase aberrations are caused by surface errors. 
lam0_ref = mp.full.lambdasMat(mp.si_ref,mp.wi_ref);
phz = angle(mp.P1.full.E(:,:,mp.wi_ref,mp.si_ref))*lam0_ref/lambda;
mp.wfs.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex) = exp(1i*phz);


%--Set the point source as the exoplanet or the star
if strcmpi(modvar.whichSource, 'exoplanet') %--Don't include tip/tilt jitter for planet wavefront since the effect is minor
    %--The planet does not move in sky angle, so the actual tip/tilt angle needs to scale inversely with wavelength.
    planetAmp = sqrt(mp.c_planet);  % Scale the E field to the correct contrast
    planetPhase = (-1)*(2*pi*(mp.x_planet*mp.P2.full.XsDL + mp.y_planet*mp.P2.full.YsDL));
    Ein = planetAmp*exp(1i*planetPhase*mp.lambda0/lambda);

elseif strcmpi(modvar.whichSource,'offaxis') %--Use for throughput calculations 
    TTphase = (-1)*(2*pi*(modvar.x_offset*mp.P2.full.XsDL + modvar.y_offset*mp.P2.full.YsDL));
    Ett = exp(1i*TTphase*mp.lambda0/lambda);
    Ein = Ett.*mp.wfs.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex); 
        
else % Default to using the starlight
    %--Include the tip/tilt in the input stellar wavefront
    if(isfield(mp,'ttx'))  % #NEWFORTIPTILT
        %--Scale by lambda/lambda0 because ttx and tty are in lambda0/D
        x_offset = mp.ttx(modvar.ttIndex)*(mp.lambda0/lambda);
        y_offset = mp.tty(modvar.ttIndex)*(mp.lambda0/lambda);

        TTphase = (-1)*(2*pi*(x_offset*mp.P2.full.XsDL + y_offset*mp.P2.full.YsDL));
        Ett = exp(1i*TTphase*mp.lambda0/lambda);
        Ein = Ett.*mp.wfs.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex);  

    else %--Backward compatible with code without tip/tilt offsets in the Jacobian
        Ein = mp.wfs.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex);  
    end
end

%--Apply a Zernike (in amplitude) at input pupil if specified
if(isfield(modvar,'zernIndex')==false)
    modvar.zernIndex = 1;
end

if(modvar.zernIndex~=1)
    indsZnoll = modvar.zernIndex; %--Just send in 1 Zernike mode
    zernMat = falco_gen_norm_zernike_maps(mp.P1.full.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    Ein = Ein.*zernMat*(2*pi/lambda)*mp.jac.Zcoef(modvar.zernIndex);
end

%%

mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.compact.NdmPad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Compute the DM surfaces for the current DM commands

if(any(mp.dm_ind==1)); DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, NdmPad); else; DM1surf = 0; end %--Pre-compute the starting DM1 surface
if(any(mp.dm_ind==2)); DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, NdmPad); else; DM2surf = 0; end %--Pre-compute the starting DM2 surface

pupil = padOrCropEven(mp.P1.compact.mask,NdmPad);
Ein = padOrCropEven(Ein,mp.compact.NdmPad);

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.compact.mask, NdmPad); else; DM1stop = 1; end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.compact.mask, NdmPad); else; DM2stop = 1; end

if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(mp.dm_ind==1)); DM1surf = gpuArray(DM1surf); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation: 2 DMs, apodizer, binary-amplitude FPM, LS, and final focal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP2 = propcustom_2FT(EP1,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 deg.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if(abs(mp.d_P2_dm1)~=0) %--E-field arriving at DM1
    Edm1 = propcustom_PTP(EP2,mp.P2.compact.dx*NdmPad,lambda,mp.d_P2_dm1);
else
    Edm1 = EP2;
end
Edm1 = DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.compact.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if(mp.d_P2_dm1 + mp.d_dm1_dm2 == 0)
    EP2eff = Edm2;
else
    EP2eff = propcustom_PTP(Edm2,mp.P2.compact.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1));
end

if(to_input) % if the user only wants to propagate to the input pupil (i.e. through primary and DMs)
    Eout = EP2eff;
else % otherwise, go through the mask and on to the next pupil.

    %--Rotate 180 degrees to propagate to pupil P3
    EP3 = propcustom_2FT(EP2eff, mp.centering);

    %--Apply apodizer mask.
    if(mp.flagApod)
        EP3 = mp.P3.compact.mask.*padOrCropEven(EP3, mp.P3.compact.Narr); 
    end

    %--MFT from SP to FPM (i.e., P3 to F3)
    EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda, ...
        mp.P2.compact.dx, mp.wfs.mask.dxi, mp.wfs.mask.Nxi, ...
        mp.wfs.mask.deta, mp.wfs.mask.Neta, mp.centering); %--E-field incident upon the FPM

    %-- Propagate from FPM to WFS cam 
    phzSupport = mp.wfs.mask.phzSupport; % support of the WFS phase dimple 
    maskDepth_m = mp.wfs.mask.depth; % mask depth in meters 
    maskAmp = mp.wfs.mask.amp; % dimple amplitude 
    WFScam_Narr = mp.wfs.cam.Narr; % Array size at WFS camera 
    WFScam_dx = mp.wfs.cam.dx;
    
    if(refwave)
        EF3 = phzSupport.*EF3inc; % Take only the part of the beam in the phase dimple
    else
        if(strcmpi(mp.wfs.mask.type,'transmissive'))
            n_mask = mp.wfs.mask.n(lambda);
            FPM = maskAmp.*exp(1i*2*pi/lambda*(n_mask-1)*maskDepth_m.*phzSupport);
        elseif(strcmpi(mp.wfs.mask.type,'reflective'))
            FPM = maskAmp.*exp(1i*4*pi/lambda*maskDepth_m.*phzSupport);
        else
            disp('mp.wfs.mask.type must be transmissive or reflective');
        end
        EF3 = (1-FPM).*EF3inc; %--Apply (1-FPM) for Babinet's principle later

        %--Use Babinet's principle at the Lyot plane.
        EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
        EP4noFPM = padOrCropEven(EP4noFPM,WFScam_Narr); %--Crop down to the size of the Lyot stop opening
    end

    %--MFT from WFS FPM to WFS camera
    EP4sub = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.wfs.mask.dxi,mp.wfs.mask.deta,WFScam_dx,WFScam_Narr,mp.centering); % Subtrahend term for Babinet's principle     

    if(refwave)
        Eout = EP4sub;
    else
        Eout = EP4noFPM-EP4sub;
    end
end

if(mp.useGPU)
    Eout = gather(Eout);
end

end % End of function