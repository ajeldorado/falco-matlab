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


function Eout = model_compact_scale(mp, lambda, Ein, normFac, flagEval, flagUseFPM)

mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.compact.NdmPad;
transOuterFPM = 1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Compute the DM surfaces for the current DM commands
if(any(mp.dm_ind==1)); DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, NdmPad); else; DM1surf = 0; end %--Pre-compute the starting DM1 surface
if(any(mp.dm_ind==2)); DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, NdmPad); else; DM2surf = 0; end %--Pre-compute the starting DM2 surface

pupil = padOrCropEven(mp.P1.compact.mask,NdmPad);
Ein = padOrCropEven(Ein,mp.compact.NdmPad);

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.compact.mask, NdmPad); else; DM1stop = ones(NdmPad); end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.compact.mask, NdmPad); else; DM2stop = ones(NdmPad); end

if mp.useGPU
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(mp.dm_ind==1)); DM1surf = gpuArray(DM1surf); end
end

%--Initialize as false if it doesn't exist
if(isfield(mp.full,'use_hlc_dm_patterns')==false)
    mp.full.use_hlc_dm_patterns = false;
end
%--Apply WFE to DMs 1 and 2
if(mp.full.use_hlc_dm_patterns || mp.flagDMwfe)
    if(any(mp.dm_ind==1));  Edm1WFE = exp(2*pi*1i/lambda.*padOrCropEven(mp.dm1.compact.wfe,NdmPad,'extrapval',0)); else; Edm1WFE = ones(NdmPad); end
    if(any(mp.dm_ind==2));  Edm2WFE = exp(2*pi*1i/lambda.*padOrCropEven(mp.dm2.compact.wfe,NdmPad,'extrapval',0)); else; Edm2WFE = ones(NdmPad); end
else
    Edm1WFE = ones(NdmPad);
    Edm2WFE = ones(NdmPad);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil .* Ein; %--E-field at pupil plane P1
EP2 = propcustom_relay(EP1,mp.Nrelay1to2,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 deg.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if abs(mp.d_P2_dm1) ~= 0; Edm1 = propcustom_PTP(EP2,mp.P2.compact.dx*NdmPad,lambda,mp.d_P2_dm1); else; Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = Edm1WFE.*DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.compact.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = Edm2WFE.*DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if mp.d_P2_dm1 + mp.d_dm1_dm2 == 0; EP2eff = Edm2; else; EP2eff = propcustom_PTP(Edm2,mp.P2.compact.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); end %--Back propagate to pupil P2

%--Rotate 180 degrees to propagate to pupil P3
EP3 = propcustom_relay(EP2eff,mp.Nrelay2to3,mp.centering);

%--Apply apodizer mask.
if mp.flagApod
    EP3 = mp.P3.compact.mask.*padOrCropEven(EP3, mp.P3.compact.Narr); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Select propagation based on coronagraph type   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if flagUseFPM
    
    switch upper(mp.coro)

        case{'HLC'}
            FPM = mp.FPM.mask; %--Complex transmission of the FPM
            transOuterFPM = FPM(1,1); %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
            scaleFac = lambda/mp.lambda0;

            %--Propagate to focal plane F3
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.compact.dx,scaleFac*mp.F3.compact.dxi,mp.F3.compact.Nxi,scaleFac*mp.F3.compact.deta,mp.F3.compact.Neta,mp.centering); %--E-field incident upon the FPM
            %--Apply (1-FPM) for Babinet's principle later
            EF3 = (transOuterFPM - FPM).*EF3inc; %- transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
            %--Use Babinet's principle at the Lyot plane.
            EP4noFPM = propcustom_relay(EP3,mp.Nrelay3to4,mp.centering); %--Propagate forward another pupil plane 
            EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.compact.Narr); %--Crop down to the size of the Lyot stop opening
            EP4noFPM = transOuterFPM*EP4noFPM; %--Apply the phase and amplitude change from the FPM's outer complex transmission.
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4sub = propcustom_mft_FtoP(EF3,mp.fl,lambda,scaleFac*mp.F3.compact.dxi,scaleFac*mp.F3.compact.deta,mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering); % Subtrahend term for Babinet's principle     
            EP4sub = propcustom_relay(EP4sub,mp.Nrelay3to4-1,mp.centering); %--Propagate forward more pupil planes if necessary.
            %--Babinet's principle at P4
            EP4 = (EP4noFPM-EP4sub); 

        otherwise
            error('%s is not an allowed value of mp.coro in model_compact_scale.', mp.coro);
    end

else % No FPM in beam path, so relay directly from P3 to P4.
    
    EP4 = propcustom_relay(EP3 ,mp.Nrelay3to4, mp.centering);
    EP4 = transOuterFPM * EP4;

    % Interpolate beam if Lyot plane has different resolution
    if mp.P4.compact.Nbeam ~= mp.P1.compact.Nbeam
        N1 = length(EP4);
        mag = mp.P4.compact.Nbeam / mp.P1.compact.Nbeam;
        N4 = ceil_even(mag*N1);
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
    end
    
    EP4 = padOrCropEven(EP4, mp.P4.compact.Narr);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Back to common propagation any coronagraph type   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Apply the Lyot stop
EP4 = mp.P4.compact.croppedMask.*EP4;

%--MFT to detector
EP4 = propcustom_relay(EP4,mp.NrelayFend,mp.centering); %--Get orientation of final image correct.
EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx, dxi,Nxi,deta,Neta, mp.centering);

%--Don't apply FPM if normalization value is being found
if normFac == 0
    Eout = EFend; %--Don't normalize if normalization value is being found
else
    Eout = EFend/sqrt(normFac); %--Apply normalization
end

if mp.useGPU
    Eout = gather(Eout);
end

end % End of entire function
