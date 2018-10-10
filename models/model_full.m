% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_full(mp,   modvar)
%--Full-knowledge optical model.
%    --> Not used by the estimator and controller.
%    --> Only used to create simulated intensity images.
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
% -DM = structure of DM settings
% -modvar = structure of model variables
%
%
% OUTPUTS:
% -Eout
%
% modvar structure fields (4):
% -sbpIndex
% -wpsbpIndex
% -whichSource
% -flagGenMat


function Eout = model_full(mp,   modvar,varargin)

% Set default values of input parameters
if(isfield(modvar,'sbpIndex'))
    normFac = mp.F4.full.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
elseif(isfield(modvar,'ebpIndex')) %--Entire bandpass index, out of mp.full.Nlam
    normFac = mp.F4.full.I00(modvar.ebpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
end

    %--Enable different arguments values by using varargin
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'normoff','unnorm','nonorm'} % Set to 0 when finding the normalization factor
        normFac = 0; 
        %fprintf('model_full: Not normalized.\n');
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
%--Select the type of coronagraph
switch mp.coro 
    case{'EHLC'} %--DMs, optional apodizer, extended FPM with metal and dielectric modulation and outer stop, and LS. Uses 1-part direct MFTs to/from FPM
        %--Complex transmission map of the FPM.
        ilam = (modvar.sbpIndex-1)*mp.Nwpsbp + modvar.wpsbpIndex;
        if( isfield(mp,'FPMcubeFull') )  %--Load it if stored
            FPM = mp.FPMcubeFull(:,:,ilam); 
        else %--Otherwise generate it
            FPM = falco_gen_EHLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'full'); %padOrCropEven( ,mp.dm9.NxiFPM);
        end
        
        Eout = model_full_EHLC(mp,   lambda, normFac, Ein, FPM); 
        
    case{'HLC','APHLC'} %--DMs, optional apodizer, FPM with optional metal and dielectric modulation, and LS. Uses Babinet's principle about FPM.
        %--Complex transmission map of the FPM.
        ilam = (modvar.sbpIndex-1)*mp.Nwpsbp + modvar.wpsbpIndex;
        if( isfield(mp,'FPMcubeFull') )  %--Load it if stored
            FPM = mp.FPMcubeFull(:,:,ilam);
        else %--Otherwise generate it
            FPM = falco_gen_HLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'full'); %padOrCropEven( ,mp.dm9.NxiFPM);
        end

        Eout = model_full_HLC(mp,   lambda, normFac, Ein, FPM);  
        
    case{'SPHLC','FHLC'} %--DMs, optional apodizer, complex/hybrid FPM with outer diaphragm, LS. Uses 2-part direct MFTs to/from FPM
        Eout = model_full_SPHLC(mp,   lambda, Ein, normFac);
        
        
    case{'LC','DMLC','APLC'} %--optional apodizer, occulting spot FPM, and LS.
        Eout = model_full_LC(mp,   lambda, Ein, normFac);
       
    case{'SPLC','FLC'} %--Optional apodizer, binary-amplitude FPM with outer diaphragm, LS
        Eout = model_full_SPLC(mp,   lambda, Ein, normFac);
            
    case{'vortex','Vortex','VC','AVC'} %--Optional apodizer, vortex FPM, LS
        Eout = model_full_VC(mp,   lambda, Ein, normFac);       
          
%     case{'SPC','APP','APC'} %--Pupil-plane mask only
%         Eout = model_full_APC(mp,   modvar);   
        
    otherwise
        error('model_full.m: Modely type\t %s\t not recognized.\n',mp.coro);        
                    
end

if(mp.useGPU)
    Eout = gather(Eout);
end

end % End of function


    
