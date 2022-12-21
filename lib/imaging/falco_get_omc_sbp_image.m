% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to:
%    set DMs
%    set wavelength band
%    get an image in the specified sub-bandpass
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
%      uses: mp.dm1.V, mp.dm2.V, mp.tb, .debug, .isProbing, .sbp_centers
% - si = index of sub-bandpass for which to take the image
%
% OUTPUTS
% - normI: Normalized intensity in the sub-bandpass
%          	(i.e. approximate raw contrast but normalized 
%           by a photometry measurement at a single offset)
%
% REVISION HISTORY
% - Modified from falco_get_gpct_sbp_image on 2019-06-26 by G. Ruane
% - Modified from falco_get_hcst_sbp_image on 2019-03-22 by G. Ruane
% - Created on 2019-03-22 by G. Ruane 
% - Copied from falco_get_dst_sbp_image.m and modified 2022-02 by D. Marx

function normI = falco_get_omc_sbp_image(mp,si)

    % convenience:
    tb = mp.tb;
    sbp_width = tb.info.sbp_width(si); %--Width of each sub-bandpass on testbed (meters)
    star_power = tb.info.star_power(si); %--star setting
    
    % already set in config: tb.sciCam.subdir = 'falco';
    if isfield(mp,'debug'), debug = mp.debug; else, debug = false; end
    NM = 1e-9;
    
    if(mp.isProbing)
        sbp_texp  = tb.info.sbp_texp_probe(si);% Exposure time for each sub-bandpass (seconds)
    else
        % TO DO: Add the capability to make the exposure time adaptive 
        % if(tb.info.adaptive_texp)
        %   query current Inorm (c = mean normI)
        %   sbp_texp = sciCam_getAdaptiveExposureTime(c,1e-8,10)
        % else
        sbp_texp  = tb.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
    end
    
    PSFpeak   = tb.info.PSFpeaks(si);% counts per second 
    
    
    %----- Send commands to the DM -----
    %disp('Sending current DM voltages to testbed') 
    % dm.V is relative to flatmap or biasMap
    if(mp.dm1.transp)
        dm1_map = mp.dm1.V'; % There's a transpose between Matlab and DM indexing
    else
        dm1_map = mp.dm1.V;
    end
    if(mp.dm2.transp)
        dm2_map = mp.dm2.V'; % There's a transpose between Matlab and DM indexing
    else
        dm2_map = mp.dm2.V;
    end

    if false,
        hfig = figure(700);
        figure_mxn(hfig, 2, 2);

        subplot(2,2,1)
        imagesc(dm1_map+tb.DM1.flatmap);
        axis image;
        colorbar;
        title('DM1 Total')
        
        subplot(2,2,2)
        imagesc(dm1_map);
        axis image;
        colorbar;
        title('DM1 - flatmap')
        drawnow;

        subplot(2,2,3)
        imagesc(dm2_map+tb.DM2.flatmap);
        axis image;
        colorbar;
        title('DM2 Total')
        
        subplot(2,2,4)
        imagesc(dm2_map);
        axis image;
        colorbar;
        title('DM2 - flatmap')
        drawnow;
    end

    % Send the commands to the DM. 
    % Note: tb.DM.flatmap contains the commands to flatten the DM. 
    %       mp.dm1.V is added to the flat commands inside DM_apply2Dmap. 
    if(tb.DM1.installed && tb.DM1.CONNECTED)
        DM_apply2Dmap(tb.DM1,dm1_map); % relative to DM.flatmap
    end
    if(tb.DM2.installed && tb.DM2.CONNECTED)
        try
            DM_apply2Dmap(tb.DM2,dm2_map); % relative to DM.flatmap
        catch ME
            error(ME);
    
            %%%% DST only:
            %             %try; cleanUpDMs(tb); end
            %             disp('Error setting DM2. Reseting electronics. Trying again.')
            %             % push the button on the DM controller (DST only)
            %             FNGR_setPos(tb,5);FNGR_setPos(tb,8);FNGR_setPos(tb,5);
            %             pause(5);
            %             setUpDMs(tb);
            %             DM_apply2Dmap(tb.DM1,dm1_map);
            %             DM_apply2Dmap(tb.DM2,dm2_map);
        end
            
    end
    
    %----- Get image from the testbed -----
    disp(['Getting image from testbed in band ',num2str(si),'. texp = ',num2str(sbp_texp)])
    
    % Set wavelength
    %disp(['Setting varia to bandpass',num2str(si)])
    lam0 = mp.sbp_centers(si);
    lam1 = lam0 - sbp_width/2;
    lam2 = lam0 + sbp_width/2;
    tb.star.lower = lam1/NM;
    tb.star.upper = lam2/NM;
    tb.star.power = star_power; 
    
    if false %debug
        disp(tb.star);
    end
    
    %     if(strcmpi(tb.info.source,'nkt'))
    %         NKT_setWvlRange(tb,lam1*1e9,lam2*1e9); % DST/gruane_DST/tb_lib/NKT/NKT_setWvlRange
    %     end
    
    % Load a dark
    dark = sciCam_loadDark(tb,sbp_texp); % DST/gruane_DST/tb_lib/scicam/sciCam_loadDark
    
    % Scale the PSF photometry by the current integration time
    PSFpeak_counts = PSFpeak*sbp_texp; 
    
    % Get normalized intensity (dark subtracted and normalized by PSFpeak)
    % sciCam_getImage returns FOV window to match falco expected image size   

    rawIm = sciCam_getImage(tb,sbp_texp);
    normI = (rawIm-dark)/PSFpeak_counts; % DST/gruane_DST/tb_lib/scicam/sciCam_getImage
    
end

