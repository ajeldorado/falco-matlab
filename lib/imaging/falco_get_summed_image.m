% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Isum = falco_get_summed_image(mp)
%
% Function to get a broadband image over the entire bandpass by summing the 
% sub-bandpass images.
%
%--INPUTS
% mp = structure of all model parameters
%
%--OUTPUTS
% Imean = 2-D stellar PSF averaged over sub-bandpasses and polarization states [normalized intensity]
%
%--REVISION HISTORY
% - Modified on 2019-05-06 by A.J. Riggs to include a different option for
% looping over a full model in PROPER.
% - Simplified on 2019-03-01 by A.J. Riggs to loop over falco_get_sbp_image.m 
%--------------------------------------------------------------------------

function [Imean,varargout] = falco_get_summed_image(mp)

    
    if(mp.flagParfor && mp.flagSim) %(mp.flagSim && mp.full.flagPROPER) %--Save a lot of time by making all PROPER full model in parallel
        %--Compute the DM surfaces outside the full model to save some time
        if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
        if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
        if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end
        %--Loop over all wavelengths and polarizations        
        ind_list = allcomb(1:mp.full.NlamUnique,1:length(mp.full.pol_conds)).';
        Nval = size(ind_list,2);
        
        %--Obtain all the images in parallel
        parfor ic=1:Nval
            if mp.flagFiber
                [Iall{ic},Ifib{ic}] = falco_get_single_sim_image(ic,ind_list,mp);  
            else
                Iall{ic} = falco_get_single_sim_image(ic,ind_list,mp);  
            end
        end

        %--Apply the spectral weights and add together
        Imean = 0;
        if mp.flagFiber;Ifibmean = 0;end
        for ic=1:Nval  
            ilam = ind_list(1,ic);
            Imean = Imean + mp.full.lambda_weights_all(ilam)/length(mp.full.pol_conds)*Iall{ic};  
            if mp.flagFiber; Ifibmean = Ifibmean + mp.full.lambda_weights_all(ilam)/length(mp.full.pol_conds)*Ifib{ic};  end
        end
        if mp.flagFiber; varargout{1} = Ifibmean; end
        
    else %--Otherwise, just loop over the function to get the sbp images. Need to do this for testbeds
        Imean = 0; % Initialize image
        if mp.flagFiber
            Ifibmean = 0;
            Ifibmean_arr = zeros(1,mp.Nsbp);
        else
            Imean_arr = zeros(1,mp.Nsbp);
        end
        for si=1:mp.Nsbp  
            if mp.flagFiber
                [Imsbp,Ifiber] = falco_get_sbp_image(mp,si);
                Ifibmean = Ifibmean + mp.sbp_weights(si)*Ifiber;
                Ifibmean_arr(si) = Ifiber;
            else
                Imsbp = falco_get_sbp_image(mp,si);
                Imean_arr(si) = mean(Imsbp(mp.Fend.corr.maskBool));
            end
            Imean = Imean +  mp.sbp_weights(si)*Imsbp;
        end
        if mp.flagFiber
            varargout{1} = Ifibmean; 
            varargout{2} = Ifibmean_arr; 
        else
            varargout{1} = []; 
            varargout{2} = Imean_arr; 
        end
    end

end %--END OF FUNCTION
