% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get a simulated image in the specified sub-bandpass.
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - si = index of sub-bandpass for which to take the image
%
% OUTPUTS
% - Isum: sub-bandpass image
%
% REVISION HISTORY
% - Modified on 2019-06-14 by A.J. Riggs to use parfor over the wavelengths
% and polarization states.
% - Created on 2019-02-06 by A.J. Riggs.

function [Isbp,varargout] = falco_get_sim_sbp_image(mp,si)

%--Compute the DM surfaces outside the full model to save lots of time
if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end

%--Number of polarization states used
Npol = length(mp.full.pol_conds);  

%--Loop over all wavelengths and polarizations        
inds_list = allcomb(1:mp.Nwpsbp,1:Npol).'; %--dimensions: [2 x mp.Nwpsbp*Npol ]
Nvals = size(inds_list,2);
    
if(mp.flagParfor) %--Save a lot of time by running full model in parallel
    parfor ic=1:Nvals
        if mp.flagFiber
            [Iall{ic},Ifib{ic}] = falco_get_single_sbp_image_WvlPol(ic,inds_list,si,mp);  
        else
            Iall{ic} = falco_get_single_sbp_image_WvlPol(ic,inds_list,si,mp);  
        end
    end %--Obtain all the images in parallel
else %--Run in serial
    for ic=Nvals:-1:1 
        if mp.flagFiber
            [Iall{ic},Ifib{ic}] = falco_get_single_sbp_image_WvlPol(ic,inds_list,si,mp);  
        else
            Iall{ic} = falco_get_single_sbp_image_WvlPol(ic,inds_list,si,mp);  
        end
    end
end

%--Apply the spectral weights and sum
Isbp = 0; 
Ifib_sbp = 0;
for ic=1:Nvals  
    Isbp = Isbp + Iall{ic};  
    if mp.flagFiber
        Ifib_sbp = Ifib_sbp + Ifib{ic};
    end
end
varargout{1} = Ifib_sbp;
end %--END OF FUNCTION

%--Function to return the weighted, normalized intensity image at a given
% wavelength in the specified sub-bandpass.
function [Iout,varargout] = falco_get_single_sbp_image_WvlPol(ic,inds_list,si,mp)

    wi = inds_list(1,ic);
    ipol = inds_list(2,ic);

    %--Get the starlight image
    modvar.sbpIndex   = si;
    modvar.wpsbpIndex = wi;
    mp.full.polaxis = mp.full.pol_conds(ipol);
    modvar.whichSource = 'star';
    if mp.flagFiber
        [Estar,Efiber] = model_full(mp, modvar);
        varargout{1} = (abs(Efiber).^2);
    else
        Estar = model_full(mp, modvar);
    end
    Iout = (abs(Estar).^2); %--Apply spectral weighting outside this function

    %--Optionally include the planet PSF
    if(mp.planetFlag)
        modvar.whichSource = 'exoplanet';
        Eplanet = model_full(mp,modvar);
        Iout = Iout + abs(Eplanet).^2; %--Apply spectral weighting outside this function
    end

    %--Apply weight within the sub-bandpass. Assume polarizations are evenly weigted.
    Iout = mp.full.lambda_weights(wi)/length(mp.full.pol_conds)*Iout;
   
end %--END OF FUNCTION