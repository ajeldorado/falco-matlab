% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass from the DST-2. 
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

function normI = falco_get_dst2_sbp_image(mp,si)

    tb = mp.tb;
    sbp_width = tb.info.sbp_width(si); %--Width of each sub-bandpass on testbed (meters)
    tb.sciCam.subdir = 'falco';
    
    if mp.isProbing
        sbp_texp  = tb.info.sbp_texp_probe(si);% Exposure time for each sub-bandpass (seconds)
    else
        sbp_texp  = tb.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
    end
    
    PSFpeak   = tb.info.PSFpeaks(si);% counts per second 
    
    
    %----- Send commands to the DM -----
    %disp('Sending current DM voltages to testbed') 
    if mp.dm1.transp 
        dm1_map = mp.dm1.V'; % Applies a transpose between Matlab and DM indexing
    else
        dm1_map = mp.dm1.V;
    end
%     if(mp.dm2.transp)
%         dm2_map = mp.dm2.V'; % There's a transpose between Matlab and DM indexing
%     else
%         dm2_map = mp.dm2.V;
%     end
    
%     figure(700)
%     subplot(1,2,1)
%     imagesc(dm1_map+tb.DM1.flatmap);
%     axis image; 
%     colorbar;
% 
%     subplot(1,2,2)
%     imagesc(dm1_map);
%     axis image; 
%     colorbar;
%     drawnow; 

    % Send the commands to the DM. 
    % Note: tb.DM.flatmap contains the commands to flatten the DM. 
    %       mp.dm1.V is added to the flat commands inside DM_apply2Dmap. 
    if tb.DM1.installed && tb.DM1.CONNECTED 
        
        try
            DM_apply2Dmap(tb.DM1,dm1_map);
        catch 
            try; cleanUpDMs(tb); end
            disp('Error setting DM1. Reseting electronics. Trying again.')
            FNGR_setPos(tb,tb.FNGR.inpos);FNGR_setPos(tb,tb.FNGR.outpos);
            pause(5);
            setUpDMs(tb);
            DM_apply2Dmap(tb.DM1,dm1_map);
            DM_apply2Dmap(tb.DM1,dm1_map);
            DM_apply2Dmap(tb.DM1,dm1_map);
        end
    end
    if tb.DM2.installed && tb.DM2.CONNECTED 
        DM_apply2Dmap(tb.DM2,dm2_map);
    end
    
    %----- Get image from the testbed -----
    disp(['Getting image from testbed in band ',num2str(si),'. texp = ',num2str(sbp_texp)])
    
    % Set wavelength
    %disp(['Setting varia to bandpass',num2str(si)])
    lam0 = mp.sbp_centers(si);
    lam1 = lam0 - sbp_width/2;
    lam2 = lam0 + sbp_width/2;
    if(strcmpi(tb.info.source,'nkt'))
        NKT_setWvlRange(tb,lam1*1e9,lam2*1e9);
    end
    
    % Load a dark
    dark = sciCam_loadDark(tb,sbp_texp);
    
    % Scale the PSF photometry by the current integration time
    PSFpeak_counts = PSFpeak*sbp_texp; 
    
    % Get normalized intensity (dark subtracted and normalized by PSFpeak)
    normI = (sciCam_getImage(tb,sbp_texp)-dark)/PSFpeak_counts; 
    
end
