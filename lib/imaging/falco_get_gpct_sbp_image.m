% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass from the GPCT. 
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
% - Modified from falco_get_hcst_sbp_image on 2019-03-22 by G. Ruane
% - Created on 2019-03-22 by G. Ruane 

function [normI,newV] = falco_get_gpct_sbp_image(mp,si)

    bench = mp.bench;
    sbp_width = bench.info.sbp_width(si); %--Width of each sub-bandpass on testbed (meters)
    sbp_texp  = bench.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
    PSFpeak   = bench.info.PSFpeaks(si);% counts per second 

    
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
    
    
    %----- Send commands to the DM -----
    %disp('Sending current DM voltages to testbed') 
    
    if(mp.dm1.run_dmsmooth)
        [newV,message] = tb_DM_dmsmooth( bench, mp.dm1.V );
    else
        newV = mp.dm1.V;
    end

    if(mp.dm1.transp)
        map = newV'; % There's a transpose between Matlab and DM indexing
    else
        map = newV;
    end

    % Send the commands to the DM. 
    % Notes: bench.DM.flatmap contains the commands to flatten the DM. 
    %        mp.dm1.V is added to the flat commands inside tb_DM_apply2Dmap. 
    tb_DM_apply2Dmap(bench,map);
    
    %----- Get image from the testbed -----
    disp(['Getting image from testbed in band ',num2str(si)])
    
    % Set wavelength
    %disp(['Setting varia to bandpass',num2str(si)])
    lam0 = mp.sbp_centers(si);
    lam1 = lam0 - sbp_width/2;
    lam2 = lam0 + sbp_width/2;
    tb_NKT_setWvlRange(bench,lam1*1e9,lam2*1e9);
    
    % Load a dark
    dark = tb_andor_loadDark(bench,sbp_texp);
    
    % Scale the PSF photometry by the current integration time
    PSFpeak_counts = PSFpeak*sbp_texp; 
    
    % Get normalized intensity (dark subtracted and normalized by PSFpeak)
    normI = (tb_andor_getImage(bench,sbp_texp)-dark)/PSFpeak_counts; 
    
end
