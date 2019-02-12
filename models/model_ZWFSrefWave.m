

function Eout = model_ZWFSrefWave(mp,   modvar,varargin)


% %--Save a lot of RAM (remove compact model data from full model inputs)
% if(any(mp.dm_ind==1)); mp.dm1 = rmfield(mp.dm1,'compact'); end

%--Set the wavelength
if(isfield(modvar,'lambda'))
    lambda = modvar.lambda;
elseif(isfield(modvar,'ebpIndex'))
    lambda = mp.full.lambdas(modvar.ebpIndex);
elseif(isfield(modvar,'sbpIndex'))
    lambda = mp.sbp_centers(modvar.sbpIndex)*mp.full.sbp_facs(modvar.wpsbpIndex);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input E-fields
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--Set the point source as the exoplanet or the star
if strcmpi(modvar.whichSource, 'exoplanet') %--Don't include tip/tilt jitter for planet wavefront since the effect is minor
    %--The planet does not move in sky angle, so the actual tip/tilt angle needs to scale inversely with wavelength.
    planetAmp = sqrt(mp.c_planet);  % Scale the E field to the correct contrast
    planetPhase = (-1)*(2*pi*(mp.x_planet*mp.P2.full.XsDL + mp.y_planet*mp.P2.full.YsDL));
    Ein = planetAmp*exp(1i*planetPhase*mp.lambda0/lambda);

elseif strcmpi(modvar.whichSource,'offaxis') %--Use for throughput calculations 
    TTphase = (-1)*(2*pi*(modvar.x_offset*mp.P2.full.XsDL + modvar.y_offset*mp.P2.full.YsDL));
    Ett = exp(1i*TTphase*mp.lambda0/lambda);
    Ein = Ett.*mp.P1.full.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex); 
        
else % Default to using the starlight
    %--Include the tip/tilt in the input stellar wavefront
    if(isfield(mp,'ttx'))  % #NEWFORTIPTILT
        %--Scale by lambda/lambda0 because ttx and tty are in lambda0/D
        x_offset = mp.ttx(modvar.ttIndex)*(mp.lambda0/lambda);
        y_offset = mp.tty(modvar.ttIndex)*(mp.lambda0/lambda);

        TTphase = (-1)*(2*pi*(x_offset*mp.P2.full.XsDL + y_offset*mp.P2.full.YsDL));
        Ett = exp(1i*TTphase*mp.lambda0/lambda);
        Ein = Ett.*mp.P1.full.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex);  

    else %--Backward compatible with code without tip/tilt offsets in the Jacobian
%         Ein = mp.Estar(:,:,modvar.wpsbpIndex,modvar.sbpIndex);  
        Ein = mp.P1.full.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex);  
    end
end


%--Apply a Zernike (in amplitude) at input pupil if specified
if(isfield(modvar,'zernIndex')==false)
    modvar.zernIndex = 1;
end

if(modvar.zernIndex~=1)
    indsZnoll = modvar.zernIndex; %--Just send in 1 Zernike mode
    zernMat = falco_gen_norm_zernike_maps(mp.P1.full.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    % figure(1); imagesc(zernMat); axis xy equal tight; colorbar; 
    Ein = Ein.*zernMat*(2*pi/lambda)*mp.jac.Zcoef(modvar.zernIndex);
end

%%

mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.compact.NdmPad;

% dxi = mp.F4.dxi;
% Nxi = mp.F4.Nxi;
% deta = mp.F4.deta;
% Neta = mp.F4.Neta; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Compute the DM surfaces for the current DM commands
% % % mp.dm1.V = 100*(eye(mp.dm1.Nact)+rot90(eye(mp.dm1.Nact),1)); %--DEBUGGING ONLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
if( abs(mp.d_P2_dm1)~=0 ); Edm1 = propcustom_PTP(EP2,mp.P2.compact.dx*NdmPad,lambda,mp.d_P2_dm1); else; Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.compact.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if( mp.d_P2_dm1 + mp.d_dm1_dm2 == 0 ); EP2eff = Edm2; else; EP2eff = propcustom_PTP(Edm2,mp.P2.compact.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); end %--Back propagate to pupil P2

%--Rotate 180 degrees to propagate to pupil P3
EP3 = propcustom_2FT(EP2eff, mp.centering);

%--Apply apodizer mask.
if(mp.flagApod)
    EP3 = mp.P3.compact.mask.*padOrCropEven(EP3, mp.P3.compact.Narr); 
end

%--MFT from SP to FPM (i.e., P3 to F3)
EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.compact.dx,mp.F3.compact.dxi,mp.F3.compact.Nxi,mp.F3.compact.deta,mp.F3.compact.Neta,mp.centering); %--E-field incident upon the FPM

FPM = mp.F3.compact.mask.phzSupport;%mp.F3.compact.mask.amp.*exp(1i*2*pi/lambda*(mp.F3.n(lambda)-1)*mp.F3.t.*mp.F3.compact.mask.phzSupport);
EF3 = FPM.*EF3inc; %--Apply (1-FPM) for Babinet's principle later


% %--Use Babinet's principle at the Lyot plane.
% EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
% EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.compact.Narr); %--Crop down to the size of the Lyot stop opening
% 

%--MFT from FPM to Lyot Plane (i.e., F3 to P4)
EP4sub = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.compact.dxi,mp.F3.compact.deta,mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering); % Subtrahend term for Babinet's principle     


Eout = EP4sub;% EP4noFPM-EP4sub;

% figure;imagesc(abs(EP4noFPM));axis image; colorbar; set(gca,'ydir','normal');
% figure;imagesc(angle(EP4noFPM));axis image; colorbar; set(gca,'ydir','normal');
% figure;imagesc(abs(EP4sub));axis image; colorbar; set(gca,'ydir','normal');
% figure;imagesc(angle(EP4sub));axis image; colorbar; set(gca,'ydir','normal');
% 

if(mp.useGPU)
    Eout = gather(Eout);
end

end % End of function


    
