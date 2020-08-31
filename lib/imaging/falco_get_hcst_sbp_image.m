% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass from the Caltech
% HCST testbed. This function will need to be replaced in order to run on a
% different testbed. Note that the number of pixels per lambda*F# is
% predetermined. 
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - si = index of sub-bandpass for which to take the image
%
% OUTPUTS
% - normI: Normalized intensity in the sub-bandpass
%          	(i.e. approximate raw contrast but normalized 
%           by a photometry measurement at a single offset)
%
% REVISION HISTORY
% - Created on 2019-02-22 by G. Ruane 

function [normI,varargout] = falco_get_hcst_sbp_image(mp,si)
    
    bench = mp.bench;
    sbp_width = bench.info.sbp_width(si); %--Width of each sub-bandpass on testbed (meters)
    sbp_texp  = bench.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
    PSFpeak   = bench.info.PSFpeaks(si);% counts per second 
    
    
    %----- Send commands to the DM -----
    disp('Sending current DM voltages to testbed') 

    
    %--Quantization of DM actuation steps based on least significant bit of the
    % DAC (digital-analog converter). In height, so called HminStep_tb
    % If HminStep_tb (minimum step in H) is defined, then quantize the DM voltages
    % This is for simulating the LSB effect (Dan's experiment) 
    if(isfield(mp.dm1,'HminStep_tb') && ~any(isnan(mp.dm1.HminStep_tb(:))))

        % If desired method is not defined, set it to the default. 
        if(~isfield(mp.dm1,'HminStepMethod')); mp.dm1.HminStepMethod = 'round'; end

        % Discretize/Quantize the DM voltages (creates dm.Vquantized)
        mp.dm1 = falco_discretize_dm_surf(mp.dm1, mp.dm1.HminStepMethod, 'tb');
        mp.dm1.V = mp.dm1.Vquantized;

    end
    
    
    map = (mp.dm1.V'); % There's a transpose between Matlab and BMC indexing
%     map = flipud(mp.dm1.V'); % There's a transpose between Matlab and BMC indexing

    if ~mp.flagFiber;map = fliplr(map);end
    % Send the commands to the DM. 
    % Notes: bench.DM.flatvec contains the commands to flatten the DM. 
    %        mp.dm1.V is added to the flat commands inside
    %        hcst_DM_apply2Dmap. 
    %        FALCO allows for searches through multiplicative gains which
    %        are handled outside of this function. 
    cmds = hcst_DM_apply2Dmap(bench,map,1);% Returns actual DM commands 
    
    if(isfield(bench.info,'source') && strcmpi(bench.info.source,'nkt'))
        %----- Get image from the testbed -----
        disp(['Getting image from testbed in band ',num2str(si)])
    
        % Set wavelength
        lam0 = mp.sbp_centers(si);
        lam1 = lam0 - sbp_width/2;
        lam2 = lam0 + sbp_width/2;
        hcst_NKT_setWvlRange(mp,bench,lam1*1e9,lam2*1e9,si);
    else
        disp('Getting image from testbed (using laser source)')
    end

    % Load the dark with the correct tint. It must exist in the dark
    % library. 
    dark = hcst_andor_loadDark(bench,[bench.info.path2darks,'dark_tint',num2str(bench.andor.tint,2),'_coadds1.fits']);
    
    % Scale the PSF photometry by the current integration time
    peakPSF = PSFpeak/mp.peakPSFtint*bench.andor.tint*mp.NDfilter_cal; 
    
    % Take image
    if mp.flagFiber && mp.flagUseCamera4EFCSMF
        hcst_andor_setSubwindow(bench,bench.andor.FocusRow,...
            bench.andor.FocusCol,bench.andor.AOIHeight);
    end

    Im = hcst_andor_getImage(bench);
    
    % Get normalized intensity (dark subtracted and normalized by peakPSF)
    normI = (Im-dark)/peakPSF; 
%     normI = fliplr((hcst_andor_getImage(bench)-dark)/peakPSF); 
        
    % Check exposure time
    if si==1 && bench.andor.tint<3 && min(normI(:))<2e-8 && ~mp.flagFiber && ~mp.efc_imageSharpening
        dh = Im(mp.Fend.corr.maskBool);
        if median(dark(:))>min(dh(:))
            disp('Exposure time is too low; adding 0.5 sec')
            hcst_andor_setExposureTime(bench,bench.andor.tint+0.5);
            Im = hcst_andor_getImage(bench);
            normI = (Im-dark)/peakPSF;
            
            % We have to change the mp.est.fudge:
            mp.est.probe.gainFudge = interp1([0.5,3],[1,0.5],bench.andor.tint); % empirically found for the 780nm laser
        end
    end
    
    %% Fiber
    if mp.flagFiber
        SMFInt0   = bench.info.SMFInt0s(si);% 
        
        if ~mp.flagUseCamera4EFCSMF
            pause(0.5)
            Vsmf = hcst_readFemtoOutput_adaptive_inV(bench,bench.Femto.averageNumReads);
            Ifiber = Vsmf/SMFInt0;
            if Ifiber<0
                while Ifiber<0
                    disp('Negative reading out of the Femto')
                    averageNumReads =  hcst_fiu_computeNumReadsNeeded(bench,1e-9);
                    Vsmf = hcst_readFemtoOutput_adaptive_inV(bench,averageNumReads);
                    Ifiber = Vsmf/SMFInt0;
                end
            elseif(Ifiber<1.3e-8 && bench.Femto.averageNumReads<1000)
                averageNumReads =  hcst_fiu_computeNumReadsNeeded(bench,Ifiber);
                Vsmf = hcst_readFemtoOutput_adaptive_inV(bench,averageNumReads);
                Ifiber = Vsmf/SMFInt0;
                while Ifiber<0
                    disp('Negative reading out of the Femto')
                    averageNumReads =  hcst_fiu_computeNumReadsNeeded(bench,1e-9);
                    Vsmf = hcst_readFemtoOutput_adaptive_inV(bench,averageNumReads);
                    Ifiber = Vsmf/SMFInt0;
                end
            end            
        else
            dark4EFCSMF = hcst_andor_loadDark(bench,[bench.info.path2darks,'dark_tint',num2str(bench.andor.tint,2),'_coadds1.fits']);
            hcst_andor_setSubwindow(bench,bench.andor.FEURow,...
                bench.andor.FEUCol,128);
            im = hcst_andor_getImage(bench)-dark4EFCSMF;
            Vsmf = hcst_fiu_aperturePhotometryOnAndor(bench,im);%max(im(:));
            Ifiber = Vsmf/(SMFInt0/mp.peakPSFtint*bench.andor.tint);
            
            % plot image to see how it looks lik
            figure(112)
            imagesc(im)
            axis image
            colorbar

        end

        


        varargout{1} = Ifiber;
    end
    
end 
