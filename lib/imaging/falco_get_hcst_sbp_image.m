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
    
    
    map = mp.dm1.V'; % There's a transpose between Matlab and BMC indexing
    
    % Send the commands to the DM. 
    % Notes: bench.DM.flatvec contains the commands to flatten the DM. 
    %        mp.dm1.V is added to the flat commands inside
    %        hcst_DM_apply2Dmap. 
    %        FALCO allows for searches through multiplicative gains which
    %        are handled outside of this function. 
    cmds = hcst_DM_apply2Dmap(bench,map,1);% Returns actual DM commands 
    
    if(isfield(bench.info,'source') && srtcmpi(bench.info.source,'nkt'))
        %----- Get image from the testbed -----
        disp(['Getting image from testbed in band ',num2str(si)])
    
        % Set wavelength
        lam0 = mp.sbp_centers(si);
        lam1 = lam0 - sbp_width/2;
        lam2 = lam0 + sbp_width/2;
        tb_NKT_setWvlRange(bench,lam1*1e9,lam2*1e9);
    else
        disp('Getting image from testbed (using laser source)')
    end

    % Load the dark with the correct tint. It must exist in the dark
    % library. 
    dark = hcst_andor_loadDark(bench,[bench.info.path2darks,'dark_tint',num2str(bench.andor.tint,2),'_coadds1.fits']);
    
    % Scale the PSF photometry by the current integration time
    peakPSF = mp.peakPSF/mp.peakPSFtint*bench.andor.tint*mp.NDfilter_cal; 
    
    % Get normalized intensity (dark subtracted and normalized by peakPSF)
    normI = (hcst_andor_getImage(bench)-dark)/peakPSF; 
    
end 
