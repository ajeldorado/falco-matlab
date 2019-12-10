% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_full_scale(mp, lambda, Ein, normFac)
% 
% Function to run the full-knowledge optical model and return the final
% E-field.
% - Not used by the estimator and controller.
% - Only used to create simulated intensity images.
%
% "fpm_scale" layout:  
% - Same as "Fourier" layout in that:
%   - For design or basic modeling only. 
%   - Instead use a full model with Fresnel propagations for more realistic simulations.
% - Different from "Fourier" layout because FPM scales with wavelength
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - lambda = wavelength in meters
% - Ein = 2D input E-field at entrance
% - normFac = intensity normalization factor 
%
%
% OUTPUTS:
% - Eout = 2-D electric field at final plane of optical layout
% - varargout{1}==Efiber = E-field at final plane when a single mode fiber
% is used
%
% REVISION HISTORY:
% --------------
% Modified on 2019-07-10 by A.J. Riggs to have the FPM scale with
% wavelength.
% Modified on 2019-04-18 by A.J. Riggs to use varargout for Efiber instead
% of having Efiber as a required output. 
% Modified on 2019-04-05 by A.J. Riggs to have the normalization be
%   computed by moving the source off-axis instead of removing the FPM.
% Modified on 2019-02-14 by G. Ruane to handle scalar vortex FPMs
% Modified on 2019-02-14 by A.J. Riggs to be the "Fourier" layout for all
% types of coronagraphs.
% Modified on 2017-10-17 by A.J. Riggs to have model_full.m be a wrapper. All the 
%  actual full models, including this one, have been moved to sub-routines for clarity.
% Modified by A.J. Riggs from hcil_simTestbed.m to model_full.m.
% Modified on 2015-02-18 by A.J. Riggs from hcil_model.m to hcil_simTestbed.m to include
%  extra errors in the model to simulate the actual testbed for fake images.

function [Eout, varargout] = model_full_scale(mp, lambda, Ein, normFac)

mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.full.NdmPad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--Change model values if the full model has a different value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isfield(mp,'full'))
    if(isfield(mp.full,'dm1'))
        if(isfield(mp.full.dm1,'xc'));  mp.dm1.xc = mp.full.dm1.xc;  end % x-center location of DM1 surface [actuator widths]
        if(isfield(mp.full.dm1,'yc'));  mp.dm1.yc = mp.full.dm1.yc;  end % y-center location of DM1 surface [actuator widths]
        if(isfield(mp.full.dm1,'V0'));  mp.dm1.V = mp.dm1.V + mp.full.dm1.V0;  end % Add some extra starting command to the voltages  [volts]
    end
    if(isfield(mp.full,'dm2'))
        if(isfield(mp.full.dm2,'xc'));  mp.dm2.xc = mp.full.dm2.xc;  end % x-center location of DM2 surface [actuator widths]
        if(isfield(mp.full.dm2,'yc'));  mp.dm2.yc = mp.full.dm2.yc;  end % y-center location of DM2 surface [actuator widths]
        if(isfield(mp.full.dm2,'V0'));  mp.dm2.V = mp.dm2.V + mp.full.dm2.V0;  end % Add some extra starting command to the voltages  [volts]
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(any(mp.dm_ind==1));  DM1surf = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,NdmPad); end
if(any(mp.dm_ind==2));  DM2surf = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,NdmPad); else; DM2surf=zeros(NdmPad); end

pupil = padOrCropEven(mp.P1.full.mask,NdmPad);
Ein = padOrCropEven(Ein,NdmPad);

if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(mp.dm_ind==1)); DM1surf = gpuArray(DM1surf); end
    if(any(mp.dm_ind==2)); DM2surf = gpuArray(DM2surf); end
end

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.full.mask, NdmPad); else; DM1stop = 1; end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.full.mask, NdmPad); else; DM2stop = 1; end

if(mp.flagDMwfe)
    if(any(mp.dm_ind==1));  Edm1WFE = exp(2*pi*1i/lambda.*padOrCropEven(mp.dm1.wfe,NdmPad,'extrapval',0)); else; Edm1WFE = ones(NdmPad); end
    if(any(mp.dm_ind==2));  Edm2WFE = exp(2*pi*1i/lambda.*padOrCropEven(mp.dm2.wfe,NdmPad,'extrapval',0)); else; Edm2WFE = ones(NdmPad); end
else
    Edm1WFE = ones(NdmPad);
    Edm2WFE = ones(NdmPad);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation: entrance pupil, 2 DMs, (optional) apodizer, vortex FPM, LS, and final focal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP2 = propcustom_relay(EP1,mp.Nrelay1to2,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 degrees mp.Nrelay1to2 times.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if( abs(mp.d_P2_dm1)~=0 ); Edm1 = propcustom_PTP(EP2,mp.P2.full.dx*NdmPad,lambda,mp.d_P2_dm1); else; Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = Edm1WFE.*DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.full.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = Edm2WFE.*DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if( mp.d_P2_dm1 + mp.d_dm1_dm2 == 0 ); EP2eff = Edm2; else; EP2eff = propcustom_PTP(Edm2,mp.P2.full.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); end %--Back propagate to pupil P2

%--Re-image to pupil P3
EP3 = propcustom_relay(EP2eff,mp.Nrelay2to3,mp.centering);

%--Apply the apodizer mask (if there is one)
if(mp.flagApod)
    EP3 = mp.P3.full.mask.*padOrCropEven(EP3, mp.P3.full.Narr); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Select propagation based on coronagraph type   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
switch upper(mp.coro)
    case{'VORTEX','VC','AVC'}
%         if(mp.flagApod==false)
%             EP3 = padOrCropEven(EP3,2^nextpow2(mp.P1.full.Narr)); %--Crop down if there isn't an apodizer mask
%         end
%         % Get FPM charge 
%         if(numel(mp.F3.VortexCharge)==1)
%             % single value indicates fully achromatic mask
%             charge = mp.F3.VortexCharge;
%         else
%             % Passing an array for mp.F3.VortexCharge with
%             % corresponding wavelengths mp.F3.VortexCharge_lambdas
%             % represents a chromatic vortex FPM
%             charge = interp1(mp.F3.VortexCharge_lambdas,mp.F3.VortexCharge,lambda,'linear','extrap');
%         end
%         EP4 = propcustom_mft_Pup2Vortex2Pup( EP3, charge, mp.P1.full.Nbeam/2, 0.3, 5, mp.useGPU );  %--MFTs
%         EP4 = padOrCropEven(EP4,mp.P4.full.Narr);

    case{'SPLC','FLC'}
%         %--MFT from SP to FPM (i.e., P3 to F3)
%         EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering); %--E-field incident upon the FPM
%         EF3 = mp.F3.full.mask.amp.*EF3inc; %--Apply FPM
%         %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
%         EP4 = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); %--E-field incident upon the Lyot stop

    case{'LC','APLC','RODDIER'}
%         %--MFT from apodizer plane to FPM (i.e., P3 to F3)
%         EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
%         % Apply (1-FPM) for Babinet's principle later
%         if(strcmpi(mp.coro,'roddier'))
%             FPM = mp.F3.full.mask.amp.*exp(1i*2*pi/lambda*(mp.F3.n(lambda)-1)*mp.F3.t.*mp.F3.full.mask.phzSupport);
%             EF3 = (1-FPM).*EF3inc; %--Apply (1-FPM) for Babinet's principle later
%         else
%             EF3 = (1-mp.F3.full.mask.amp).*EF3inc;
%         end
%         % Use Babinet's principle at the Lyot plane. This is the term without the FPM.
%         EP4noFPM = propcustom_relay(EP3,mp.Nrelay3to4,mp.centering); %--Propagate forward another pupil plane 
%         %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
%         EP4subtrahend = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
%         %--Babinet's principle at P4
%         EP4 = padOrCropEven(EP4noFPM,mp.P4.full.Narr) - EP4subtrahend;

    case{'HLC'}
        FPM = mp.FPM.mask; %--Complex transmission of the FPM
        transOuterFPM = FPM(1,1); %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
        scaleFac = lambda/mp.lambda0;

        %--MFT from apodizer plane to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,scaleFac*mp.F3.full.dxi,mp.F3.full.Nxi,scaleFac*mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
        % Apply (1-FPM) for Babinet's principle later
        EF3 = (transOuterFPM-mp.FPM.mask).*EF3inc; %- transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
        % Use Babinet's principle at the Lyot plane.
        EP4noFPM = propcustom_relay(EP3,mp.Nrelay3to4,mp.centering); %--Propagate forward another pupil plane 
        EP4noFPM = transOuterFPM*padOrCropEven(EP4noFPM,mp.P4.full.Narr); %--Apply the phase and amplitude change from the FPM's outer complex transmission.
        %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
        EP4subtra = propcustom_mft_FtoP(EF3,mp.fl,lambda,scaleFac*mp.F3.full.dxi,scaleFac*mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
        %--Babinet's principle at P4
        EP4 = EP4noFPM-EP4subtra;


%         %--From compact model:
%         FPM = mp.FPM.mask; %--Complex transmission of the FPM
%         transOuterFPM = FPM(1,1); %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
%         scaleFac = lambda/mp.lambda0;
%         
%         %--Propagate to focal plane F3
%         EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.compact.dx,scaleFac*mp.F3.compact.dxi,mp.F3.compact.Nxi,scaleFac*mp.F3.compact.deta,mp.F3.compact.Neta,mp.centering); %--E-field incident upon the FPM
%         %--Apply (1-FPM) for Babinet's principle later
%         EF3 = (transOuterFPM - FPM).*EF3inc; %- transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
%         %--Use Babinet's principle at the Lyot plane.
%         EP4noFPM = propcustom_relay(EP3,mp.Nrelay3to4,mp.centering); %--Propagate forward another pupil plane 
%         EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.compact.Narr); %--Crop down to the size of the Lyot stop opening
%         EP4noFPM = transOuterFPM*EP4noFPM; %--Apply the phase and amplitude change from the FPM's outer complex transmission.
%         %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
%         EP4sub = propcustom_mft_FtoP(EF3,mp.fl,lambda,scaleFac*mp.F3.compact.dxi,scaleFac*mp.F3.compact.deta,mp.P4.compact.dx,mp.P4.compact.Narr,mp.centering); % Subtrahend term for Babinet's principle     
%         EP4sub = propcustom_relay(EP4sub,mp.Nrelay3to4-1,mp.centering); %--Propagate forward more pupil planes if necessary.
%         %--Babinet's principle at P4
%         EP4 = (EP4noFPM-EP4sub); 

    otherwise
        error('model_full_Fourier.m: Modely type\t %s\t not recognized.\n',mp.coro);
end

%--Remove the FPM completely if normalization value is being found
if(normFac==0)
    switch upper(mp.coro)
        case{'VORTEX','VC','AVC'}
            EP4 = propcustom_relay(EP3,mp.Nrelay3to4, mp.centering);
            EP4 = padOrCropEven(EP4,mp.P4.full.Narr);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Back to common propagation any coronagraph type   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Apply the Lyot stop
EP4 = mp.P4.full.croppedMask.*EP4; %padOrCropEven(EP4,mp.P4.full.Narr);

%--MFT from Lyot Stop to final focal plane (i.e., P4 to Fend)
EP4 = propcustom_relay(EP4,mp.NrelayFend,mp.centering); %--Rotate the final image 180 degrees if necessary
EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.full.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering);

%--Don't apply FPM if normalization value is being found
if(normFac==0)
    Eout = EFend;  %--Don't normalize if normalization value is being found
else
    Eout = EFend/sqrt(normFac); %--Apply normalization
end

if(mp.useGPU); Eout = gather(Eout); end

if(mp.flagFiber)
    Efiber = cell(mp.Fend.Nlens,1);
    
    for nlens = 1:mp.Fend.Nlens
        EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering,'xfc',mp.Fend.x_lenslet_phys(nlens),'yfc',mp.Fend.y_lenslet_phys(nlens));
        Elenslet = EFend.*mp.Fend.lenslet.mask;
        EF5 = propcustom_mft_PtoF(Elenslet,mp.lensletFL,lambda,mp.Fend.dxi,mp.F5.dxi,mp.F5.Nxi,mp.F5.deta,mp.F5.Neta,mp.centering);
        Efiber{nlens} = mp.F5.fiberMode(:).*sum(sum(mp.F5.fiberMode.*conj(EF5)));
    end

    Efiber = permute(reshape(cell2mat(Efiber)', mp.F5.Nxi, mp.F5.Neta, mp.Fend.Nlens), [2,1,3]);
    varargout{1} = Efiber;
end


end %--END OF FUNCTION