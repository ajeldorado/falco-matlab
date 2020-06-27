% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass from the DST. 
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
% - Modified from falco_get_gpct_sbp_image on 2019-06-26 by G. Ruane
% - Modified from falco_get_hcst_sbp_image on 2019-03-22 by G. Ruane
% - Created on 2019-03-22 by G. Ruane 

function normI = falco_get_dst_sbp_image(mp,si)

    tb = mp.tb;
    sbp_width = tb.info.sbp_width(si); %--Width of each sub-bandpass on testbed (meters)
    
    if(mp.probing)
        sbp_texp  = tb.info.sbp_texp_probe(si);% Exposure time for each sub-bandpass (seconds)
    else
        sbp_texp  = tb.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
    end
    
    PSFpeak   = tb.info.PSFpeaks(si);% counts per second 
    
    
    %----- Send commands to the DM -----
    %disp('Sending current DM voltages to testbed') 
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

    % Send the commands to the DM. 
    % Note: tb.DM.flatmap contains the commands to flatten the DM. 
    %       mp.dm1.V is added to the flat commands inside DM_apply2Dmap. 
    DM_apply2Dmap(tb.DM1,dm1_map);
    DM_apply2Dmap(tb.DM2,dm2_map);
    
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
