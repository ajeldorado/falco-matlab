% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_full(mp,modvar,varargin)
%--Full-knowledge optical model.
%    --> Not used by the estimator and controller.
%    --> Used only to create simulated intensity images.
%
% REVISION HISTORY:
% --------------
% Modified on 2017-10-17 by A.J. Riggs to have model_full.m be a wrapper. All the 
%  actual full models have been moved to sub-routines for clarity.
% model_full.m - Modified from hcil_simTestbed.m
% hcil_simTestbed.m - 18 Feb 2015: Modified from hcil_model.m. Includes
%  extra errors in the model to simulate the actual testbed for fake images.
%
% ---------------
% INPUTS:
% -mp = structure of model parameters
% -modvar = structure of model variables
%
% OUTPUTS:
% -Eout
%
% modvar structure fields (4):
% -sbpIndex
% -wpsbpIndex
% -whichSource
% -flagGenMat

function [Eout, Efiber] = model_full(mp,modvar,varargin)

% Set default values of input parameters
if(isfield(modvar,'sbpIndex'))
    normFac = mp.Fend.full.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
end

%--Enable different arguments values by using varargin
icav = 0; % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'normoff','unnorm','nonorm'} % Set to 0 when finding the normalization factor
        normFac = 0; 
      otherwise
        error('model_full: Unknown keyword: %s\n', varargin{icav});
    end
end

%--Save a lot of RAM (remove compact model data from full model inputs)
if(any(mp.dm_ind==1)); mp.dm1 = rmfield(mp.dm1,'compact'); end
if(any(mp.dm_ind==2)); mp.dm2 = rmfield(mp.dm2,'compact'); end
if(any(mp.dm_ind==8)); mp.dm8 = rmfield(mp.dm8,'compact'); end
if(any(mp.dm_ind==9)); mp.dm9 = rmfield(mp.dm9,'compact'); end

%--Set the wavelength
if(isfield(modvar,'lambda'))
    lambda = modvar.lambda;
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
    Ein = planetAmp*exp(1i*planetPhase*mp.lambda0/lambda).*mp.P1.full.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex);

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
        Ein = mp.P1.full.E(:,:,modvar.wpsbpIndex,modvar.sbpIndex);  
    end
end

%--Shift the source off-axis to compute the intensity normalization value.
%  This replaces the previous way of taking the FPM out in the optical model.
if(normFac==0)
    source_x_offset = mp.source_x_offset_norm; %--source offset in lambda0/D for normalization
    source_y_offset = mp.source_y_offset_norm; %--source offset in lambda0/D for normalization
    TTphase = (-1)*(2*pi*(source_x_offset*mp.P2.full.XsDL + source_y_offset*mp.P2.full.YsDL));
    Ett = exp(1i*TTphase*mp.lambda0/lambda);
    Ein = Ett.*mp.P1.full.E(:,:,modvar.sbpIndex); 
end

%--Apply a Zernike (in amplitude) at input pupil if specified
if(isfield(modvar,'zernIndex')==false)
    modvar.zernIndex = 1;
end

if(modvar.zernIndex~=1)
    indsZnoll = modvar.zernIndex; %--Just send in 1 Zernike mode
    zernMat = falco_gen_norm_zernike_maps(mp.P1.full.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    zernMat = padOrCropEven(zernMat,mp.P1.full.Narr);
    Ein = Ein.*zernMat*(2*pi/lambda)*mp.jac.Zcoef(modvar.zernIndex);
end

%% Pre-compute the FPM first for HLC as mp.FPM.mask
switch lower(mp.layout)
    case{'fourier'}
        
        switch upper(mp.coro) 
            case{'EHLC'} %--DMs, optional apodizer, extended FPM with metal and dielectric modulation and outer stop, and LS. Uses 1-part direct MFTs to/from FPM
                %--Complex transmission map of the FPM.
                ilam = (modvar.sbpIndex-1)*mp.Nwpsbp + modvar.wpsbpIndex;
                if(isfield(mp,'FPMcubeFull')) %--Load it if stored
                    mp.FPM.mask = mp.FPMcubeFull(:,:,ilam); 
                    mp.FPM.mask = falco_gen_EHLC_FPM_complex_trans_mat(mp,modvar.sbpIndex,modvar.wpsbpIndex,'full');
                end

            case{'HLC'} %--DMs, optional apodizer, FPM with optional metal and dielectric modulation, and LS. Uses Babinet's principle about FPM.
                %--Complex transmission map of the FPM.
                ilam = (modvar.sbpIndex-1)*mp.Nwpsbp + modvar.wpsbpIndex;
                if(isfield(mp,'FPMcubeFull')) %--Load it if stored
                    mp.FPM.mask = mp.FPMcubeFull(:,:,ilam);
                else %--Otherwise generate it
                    mp.FPM.mask = falco_gen_HLC_FPM_complex_trans_mat(mp,modvar.sbpIndex,modvar.wpsbpIndex,'full');
                end
        end
end
%% Select which optical layout's full model to use.
switch lower(mp.layout)
    case{'fourier'}
        [Eout, Efiber] = model_full_Fourier(mp, lambda, Ein, normFac);
        
    case{'wfirst_phaseb_simple'} %--Use compact model as the full model.
                
        optval = mp.full;
        optval.use_dm1 = true;
        optval.use_dm2 = true;
        optval.dm1_m = mp.dm1.V.*mp.dm1.VtoH; %--DM1 commands in meters
        optval.dm2_m = mp.dm2.V.*mp.dm2.VtoH; %--DM2 commands in meters
        if(normFac==0)
            optval.source_x_offset = -mp.source_x_offset_norm;
            optval.source_y_offset = -mp.source_y_offset_norm;
        end

        Eout = prop_run('wfirst_phaseb_v2b_compact', lambda*1e6, mp.Fend.Nxi, 'quiet', 'passvalue',optval ); %--wavelength needs to be in microns instead of meters for PROPER
        Eout = circshift(rot90(Eout,2),[1,1]);	%   rotate to same orientation as FALCO
        if(normFac~=0)
            Eout = Eout/sqrt(normFac);
        end

    case{'wfirst_phaseb_proper'} %--Use the true full model as the full model

        optval = mp.full;
        optval.use_dm1 = true;
        optval.use_dm2 = true;
        optval.dm1_m = mp.dm1.V.*mp.dm1.VtoH + mp.full.dm1.flatmap; %--DM1 commands in meters
        optval.dm2_m = mp.dm2.V.*mp.dm2.VtoH + mp.full.dm2.flatmap; %--DM2 commands in meters
        if(normFac==0)
            optval.source_x_offset = -mp.source_x_offset_norm;
            optval.source_y_offset = -mp.source_y_offset_norm;
        end

        Eout = prop_run('wfirst_phaseb_v2b', lambda*1e6, mp.Fend.Nxi, 'quiet', 'passvalue',optval ); %--wavelength needs to be in microns instead of meters for PROPER
        Eout = circshift(rot90(Eout,2),[1,1]);	%   rotate to same orientation as FALCO
        if(normFac~=0)
            Eout = Eout/sqrt(normFac);
        end
end

%% Undo GPU variables if they exist
if(mp.useGPU)
    Eout = gather(Eout);
end

end % End of function