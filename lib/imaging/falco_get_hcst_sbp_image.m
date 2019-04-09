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


function normI = falco_get_hcst_sbp_image(mp,si)
    
    %----- Send commands to the DM -----
    disp('Sending current DM voltages to testbed') 
    
    map = mp.dm1.V'; % There's a transpose between Matlab and BMC indexing
    
    % Send the commands to the DM. 
    % Notes: bench.DM.flatvec contains the commands to flatten the DM. 
    %        mp.dm1.V is added to the flat commands inside
    %        hcst_DM_apply2Dmap. 
    %        FALCO allows for searches through multiplicative gains which
    %        are handled outside of this function. 
    cmds = hcst_DM_apply2Dmap(bench,map,1);% Returns actual DM commands 
    
    %----- Get image from the testbed -----
    disp(['Getting image from testbed in band',num2str(si)])
    
    % disp('Setting varia to bandpass',num2str(si)])
    %   NOT IMPLEMENTED YET!!! ONLY WILL WORK FOR MONOCHROMATIC EFC
    
    % Load the dark with the correct tint. It must exist in the dark
    % library. 
    dark = hcst_andor_loadDark(bench,[bench.info.path2darks,'dark_tint',num2str(bench.andor.tint,2),'_coadds1.fits']);
    
    % Scale the PSF photometry by the current integration time
    peakPSF = mp.peakPSF/mp.peakPSFtint*bench.andor.tint; 
    
    % Get normalized intensity (dark subtracted and normalized by peakPSF)
    normI = (hcst_andor_getImage(bench)-dark)/peakPSF; 
    
end 

