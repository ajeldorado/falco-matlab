% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_full_Fourier(mp, lambda, Ein, normFac)
% 
% Function to run the full-knowledge optical model and return the final
% E-field.
% - Not used by the estimator and controller.
% - Only used to create simulated intensity images.
%
% "Fourier" layout:  
% - For design or basic modeling only. 
% - Instead use a full model with Fresnel propagations for more realistic simulations.
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
%
%
% OUTPUTS:
% - Eout = 2-D electric field at final plane of optical layout
%
% REVISION HISTORY:
% --------------
% Modified on 2019-02-14 by G. Ruane to handle scalar vortex FPMs
% Modified on 2019-02-14 by A.J. Riggs to be the "Fourier" layout for all
% types of coronagraphs.
% Modified on 2017-10-17 by A.J. Riggs to have model_full.m be a wrapper. All the 
%  actual full models, including this one, have been moved to sub-routines for clarity.
% Modified by A.J. Riggs from hcil_simTestbed.m to model_full.m.
% Modified on 2015-02-18 by A.J. Riggs from hcil_model.m to hcil_simTestbed.m to inclue
%  extra errors in the model to simulate the actual testbed for fake images.


function Eout = model_full_Fourier(mp, lambda, Ein, normFac)


mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.full.NdmPad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(any(mp.dm_ind==1))
    if( isfield(mp.dm1,'surfM') );   DM1surf = padOrCropEven(mp.dm1.surfM, NdmPad);
    else                             DM1surf = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,NdmPad); end
else
    DM1surf = 0;
end
if(any(mp.dm_ind==2))
    if( isfield(mp.dm2,'surfM') );   DM2surf = padOrCropEven(mp.dm2.surfM, NdmPad);
    else                             DM2surf = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,NdmPad); end
else
    DM2surf = 0;
end


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
EP2 = propcustom_2FT(EP1,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 deg.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if( abs(mp.d_P2_dm1)~=0 ); Edm1 = propcustom_PTP(EP2,mp.P2.full.dx*NdmPad,lambda,mp.d_P2_dm1); else; Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = Edm1WFE.*DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.full.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = Edm2WFE.*DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if( mp.d_P2_dm1 + mp.d_dm1_dm2 == 0 ); EP2eff = Edm2; else; EP2eff = propcustom_PTP(Edm2,mp.P2.full.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); end %--Back propagate to pupil P2

%--Rotate 180 degrees to propagate to pupil P3
EP3 = propcustom_2FT(EP2eff, mp.centering);

%--Apply the apodizer mask (if there is one)
if(mp.flagApod)
    EP3 = mp.P3.full.mask.*padOrCropEven(EP3, mp.P3.full.Narr); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Select propagation based on coronagraph type   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Don't apply FPM if normalization value is being found
if(normFac==0)
    
    switch upper(mp.coro)
        case{'VORTEX','VC','AVC'}
            EP4 = propcustom_2FT(EP3, mp.centering);
            EP4 = padOrCropEven(EP4,mp.P4.full.Narr);
            
        case{'SPLC','FLC'}
            %--MFT from SP to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering); %--E-field incident upon the FPM
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4 = propcustom_mft_FtoP(EF3inc,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); %--E-field incident upon the Lyot stop 
        
        case{'LC','APLC','RODDIER'}
            EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Re-image forward (no FPM in between pupil planes)
            EP4 = mp.P4.full.croppedMask.*padOrCropEven(EP4noFPM,mp.P4.full.Narr);
           
        case{'HLC'}
            %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
            t_Ti_base = 0;
            t_Ni_vec = 0;
            t_PMGI_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
            pol = 2;
            [tCoef, ~] = falco_thin_film_material_def(lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, lambda*mp.FPM.d0fac, pol);
            transOuterFPM = tCoef;

            %--Do NOT apply FPM if normalization value is being found
            EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
            EP4 = transOuterFPM*padOrCropEven(EP4noFPM,mp.P4.full.Narr); %--Apply the phase and amplitude change from the FPM's outer complex transmission.
            
        case{'EHLC'}
            %--Complex transmission of the points outside the inner part of the FPM (just fused silica with optional dielectric and no metal).
            t_Ti_base = 0;
            t_Ni_vec = 0;
            t_PMGI_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
            pol = 2;
            [tCoef, ~] = falco_thin_film_material_def(lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, lambda*mp.FPM.d0fac, pol);
            transOuterFPM = tCoef;
            
            %--MFT from apodizer plane to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
            %--Don't apply FPM if normalization value is being found. Just
            % apply the complex transmission of the material
            EF3 = transOuterFPM.*EF3inc;
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4 = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     

        case{'FOHLC'}
            %--Do NOT apply FPM if normalization value is being found
            %--NOTE: transOuterFPM is divided out to unity in the FOHLC model
            EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
            EP4 = padOrCropEven(EP4noFPM,mp.P4.full.Narr);
            
        otherwise
            error('model_full_Fourier.m: Model type\t %s\t not recognized.\n',mp.coro);
    end
    
else
    
    switch upper(mp.coro)
        case{'VORTEX','VC','AVC'}
            if(mp.flagApod==false)
                EP3 = padOrCropEven(EP3,2^nextpow2(mp.P1.full.Narr)); %--Crop down if there isn't an apodizer mask
            end
            % Get FPM charge 
            if(numel(mp.F3.VortexCharge)==1)
                % single value indicates fully achromatic mask
                charge = mp.F3.VortexCharge;
            else
                % Passing an array for mp.F3.VortexCharge with
                % corresponding wavelengths mp.F3.VortexCharge_lambdas
                % represents a chromatic vortex FPM
                charge = interp1(mp.F3.VortexCharge_lambdas,mp.F3.VortexCharge,lambda,'linear','extrap');
            end
            EP4 = propcustom_mft_Pup2Vortex2Pup( EP3, charge, mp.P1.full.Nbeam/2, 0.3, 5, mp.useGPU );  %--MFTs
            EP4 = padOrCropEven(EP4,mp.P4.full.Narr);
            
        case{'SPLC','FLC'}
            %--MFT from SP to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering); %--E-field incident upon the FPM
            EF3 = mp.F3.full.mask.amp.*EF3inc; %--Apply FPM
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4 = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); %--E-field incident upon the Lyot stop
            
        case{'LC','APLC','RODDIER'}
            %--MFT from apodizer plane to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
            % Apply (1-FPM) for Babinet's principle later
            if(strcmpi(mp.coro,'roddier'))
                FPM = mp.F3.full.mask.amp.*exp(1i*2*pi/lambda*(mp.F3.n(lambda)-1)*mp.F3.t.*mp.F3.full.mask.phzSupport);
                EF3 = (1-FPM).*EF3inc; %--Apply (1-FPM) for Babinet's principle later
            else
                EF3 = (1-mp.F3.full.mask.amp).*EF3inc;
            end
            % Use Babinet's principle at the Lyot plane. This is the term without the FPM.
            EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4subtrahend = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
            %--Babinet's principle at P4
            EP4 = padOrCropEven(EP4noFPM,mp.P4.full.Narr) - EP4subtrahend;
            
        case{'HLC'}
            %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
            t_Ti_base = 0;
            t_Ni_vec = 0;
            t_PMGI_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
            pol = 2;
            [tCoef, ~] = falco_thin_film_material_def(lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, lambda*mp.FPM.d0fac, pol);
            transOuterFPM = tCoef;
            
            %--MFT from apodizer plane to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
            % Apply (1-FPM) for Babinet's principle later
            EF3 = (transOuterFPM-mp.FPM.mask).*EF3inc; %- transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
            % Use Babinet's principle at the Lyot plane.
            EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
            EP4noFPM = transOuterFPM*padOrCropEven(EP4noFPM,mp.P4.full.Narr); %--Apply the phase and amplitude change from the FPM's outer complex transmission.
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4subtra = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
            %--Babinet's principle at P4
            EP4 = EP4noFPM-EP4subtra;
            
        case{'EHLC'}
            %--Complex transmission of the points outside the inner part of the FPM (just fused silica with optional dielectric and no metal).
            t_Ti_base = 0;
            t_Ni_vec = 0;
            t_PMGI_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
            pol = 2;
            [tCoef, ~] = falco_thin_film_material_def(lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, lambda*mp.FPM.d0fac, pol);
            transOuterFPM = tCoef;
            
            %--MFT from apodizer plane to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
            EF3 = mp.FPM.mask.*EF3inc; %--Apply FPM
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4 = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     

        case{'FOHLC'}
            %--FPM representation (idealized as amplitude and phase)
            DM8amp = falco_gen_HLC_FPM_amplitude_from_cube(mp.dm8,'full');
            DM8ampPad = padOrCropEven( DM8amp,mp.full.Nfpm,'extrapval',1);
            DM9surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm9,'full');
            DM9surfPad = padOrCropEven( DM9surf,mp.full.Nfpm);
            transOuterFPM = 1; %--Because the complex transmission far away is divided out.
            FPM = DM8ampPad.*exp(2*pi*1i/lambda*DM9surfPad);

            %--MFT from apodizer plane to FPM (i.e., P3 to F3)
            EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
            % Apply (1-FPM) for Babinet's principle later
            EF3 = (1 - FPM/transOuterFPM).*EF3inc; %- transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
            % Use Babinet's principle at the Lyot plane.
            EP4noFPM = propcustom_2FT(EP3,mp.centering); %--Propagate forward another pupil plane 
            EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.full.Narr);
            %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
            EP4subtra = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
            %--Babinet's principle at P4
            EP4 = EP4noFPM - EP4subtra;
            
        otherwise
            error('model_full_Fourier.m: Modely type\t %s\t not recognized.\n',mp.coro);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Back to common propagation any coronagraph type   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Apply the Lyot stop
EP40 = EP4;
EP4 = mp.P4.full.croppedMask.*EP4; %padOrCropEven(EP4,mp.P4.full.Narr);


%--MFT from Lyot Stop to final focal plane (i.e., P4 to Fend.
EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.full.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta);

%--Don't apply FPM if normalization value is being found
if(normFac==0)
    Eout = EFend;  %--Don't normalize if normalization value is being found
else
    Eout = EFend/sqrt(normFac); %--Apply normalization
end

if(mp.useGPU); Eout = gather(Eout); end

if(isfield(mp,'flagElyot'))
    Eout = EP40;
end

end %--END OF FUNCTION


    
