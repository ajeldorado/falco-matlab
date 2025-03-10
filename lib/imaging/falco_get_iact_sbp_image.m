% Copyright 2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass from the IACT. 
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
% - Modified from falco_get_dst_sbp_image on 2021-04-07 by G. Ruane
% - Modified from falco_get_gpct_sbp_image on 2019-06-26 by G. Ruane
% - Modified from falco_get_hcst_sbp_image on 2019-03-22 by G. Ruane
% - Created on 2019-03-22 by G. Ruane 

function normI = falco_get_iact_sbp_image(mp,si)

    tb = mp.tb;
    sbp_width = tb.info.sbp_width(si); %--Width of each sub-bandpass on testbed (meters)
    tb.sciCam.subdir = 'falco';
    
    if mp.isProbing  % use these values if you're probing 
        sbp_texp  = tb.info.sbp_texp_probe(si);% Exposure time for each sub-bandpass (seconds)
        numCoadds = tb.info.sbp_numCoadds_probe(si);% Number of coadds to use for each sub-bandpass
    else % use these values for dark hole images 
        sbp_texp  = tb.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
        numCoadds = tb.info.sbp_numCoadds(si);% Number of coadds to use for each sub-bandpass
    end    
    tb.sciCam.numCoadds = numCoadds; % This sets the camera to use coadds
    
    % PSF photometry 
    PSFpeak   = tb.info.PSFpeaks(si);% counts per second 
    
    
	%----- Send commands to the DM -----

    % Note: tb.DM.flatmap contains the commands to flatten the DM. 
    %       mp.dm1.V is added to the flat commands inside DM_apply2Dmap. 
    if tb.DM.installed && tb.DM.CONNECTED 
        DM_apply2Dmap(tb.DM,mp.dm1.V);
    end

    
	%----- Get image from the testbed -----
    disp(['Getting image from testbed in band ',num2str(si),'. texp = ',num2str(sbp_texp)])
%     disp(['Getting image from testbed. texp=',num2str(sbp_texp),'s. numCoadds=',num2str(numCoadds)])
    
    % Set wavelength
    %disp(['Setting varia to bandpass',num2str(si)])
    lam0 = mp.sbp_centers(si);
    lam1 = lam0 - sbp_width/2;
    lam2 = lam0 + sbp_width/2;
    if strcmpi(tb.info.source,'nkt') && mp.Nsbp>1
        NKT_setWvlRange(tb,lam1*1e9,lam2*1e9);
    else
        disp('Monochromatic mode, laser not queried')
    end
    
    % Load a dark
    dark = sciCam_loadDark(tb,sbp_texp);
    
    % Scale the PSF photometry by the current integration time
    %PSFpeak_counts = PSFpeak*sbp_texp;
    PSFpeak_counts = PSFpeak*sbp_texp*numCoadds;
    
    % Get normalized intensity (dark subtracted and normalized by PSFpeak)
    normI = (sciCam_getImage(tb,sbp_texp)-dark)/PSFpeak_counts; 
    
end
