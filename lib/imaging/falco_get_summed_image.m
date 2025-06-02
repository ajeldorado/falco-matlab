% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% summedImage = falco_get_summed_image(mp)
%
% Get a broadband image over the entire bandpass by summing subband images.
%
% INPUTS
% ------
% mp : structure of all model parameters
%
% OUTPUTS
% -------
% summedImage : PSF averaged over all subbands and polarization states [normalized intensity]

function [summedImage,varargout] = falco_get_summed_image(mp)



    if(mp.flagParfor && mp.flagSim) %--Save a lot of time by making all calls to the PROPER full model in parallel

        %--Compute the DM surfaces outside the full model to save some time
        if any(mp.dm_ind == 1); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
        if any(mp.dm_ind == 2); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
        if any(mp.dm_ind == 9); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end
        
        %--Loop over all wavelengths, polarizations, and stars       
        indexCombos = allcomb(1:mp.full.NlamUnique, 1:length(mp.full.pol_conds), 1:mp.star.count).';
        Ncombos = size(indexCombos, 2);

        %--Remove testbed objects in parfor loops
        if isfield(mp, 'tb'); mp = rmfield(mp, 'tb'); end

        %--Obtain all the images in parallel
        parfor iCombo = 1:Ncombos
            if mp.flagFiber
                [Iall{iCombo},Iallfiber{iCombo}] = falco_get_single_sim_image(iCombo, indexCombos, mp);
            else
                Iall{iCombo} = falco_get_single_sim_image(iCombo, indexCombos, mp);
            end
        end
   

        %--Apply the spectral and stellar weights, and then sum
        summedImage = 0;
        summedIfiber = 0;
        for iCombo = 1:Ncombos
            iLambda = indexCombos(1, iCombo);
            iStar = indexCombos(3, iCombo);
            summedImage = summedImage + mp.star.weights(iStar)* (mp.full.lambda_weights_all(iLambda)/length(mp.full.pol_conds)) * Iall{iCombo};
            if mp.flagFiber
                summedIfiber = summedIfiber + mp.sbp_weights(iStar)* (mp.full.lambda_weights_all(iLambda)/length(mp.full.pol_conds)) * Iallfiber{iCombo};
            end
        end

    else
        
        summedImage = 0;
        summedIfiber = 0;
        for iSubband = 1:mp.Nsbp    
            if mp.flagFiber
                [sbpim,sbIfiber] = falco_get_sbp_image(mp, iSubband);
                summedIfiber = summedIfiber + mp.sbp_weights(iSubband)*sbIfiber;
                varargout{1} = summedIfiber;
            else
                sbpim = falco_get_sbp_image(mp, iSubband);
            end
            summedImage = summedImage +  mp.sbp_weights(iSubband)*sbpim;
        end
    end

end
