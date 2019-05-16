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
% Ibandavg = band-averaged image in units of normalized intensity
%
%--REVISION HISTORY
% - Modified on 2019-05-03 by A.J. Riggs from falco_get_summed_image.m to
% falco_get_summed_image_TipTiltPol.m to include tip/tilt (stellar diameter
% and pointing jitter) and differential polarization aberrations.
% - Simplified on 2019-03-01 by A.J. Riggs to loop over falco_get_sbp_image.m 
%--------------------------------------------------------------------------

function Imean = falco_get_summed_image_TipTiltPol(mp)

    %--Compute the DM surfaces outside the full model to save some time
    if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
    if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
    if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end

    Imean = 0; % Initialize image
    
    %--Number of polarization states used
    mp.full.dummy = 1; %--Initialize if this doesn't exist
    if(isfield(mp.full,'pol_conds'))  
        Npol = length(mp.full.pol_conds);  
    else
        Npol = 1;
    end
    
    %--Generate the tip/tilt offsets and their normalized weights
    [mp.full.xsTT,mp.full.ysTT,mp.full.wsTT] = falco_gen_RMS_TipTilt(mp.full.TTrms,mp.full.Dstar,mp.full.Dtel,mp.lambda0,'Nacross',mp.full.TipTiltNacross);
    Ntt = length(mp.full.wsTT);
    fprintf('\n%d tip-tilt offset points used.\n',Ntt);
    %--Iterate over all combinations of sub-bandpass, wavelength, tip/tilt offset, and polarization state.
    
    
%     InormCube = zeros(Ncombos,1);

    
    if(mp.flagParfor && mp.flagSim) %(mp.flagSim && mp.full.flagPROPER) %--Save a lot of time by making all PROPER full model in parallel
        %--Loop over all wavelengths, tip/tilt offsets, and polarizations        
        inds_list = allcomb(1:mp.full.NlamUnique,1:Ntt,1:Npol).'; %--dimensions: [3 x mp.full.NlamUnique*Ntt*Npol ]
        %inds_list = allcomb(1:mp.Nsbp,1:mp.Nwpsbp,1:Ntt,1:Npol).'; %--dimensions: [4 x mp.Nsbp*Nwpsbp*Ntt*Npol ]
        Nvals = size(inds_list,2);
        
        %--Obtain all the images in parallel
        parfor ic=1:Nvals
            Iall{ic} = falco_get_single_sim_image_TipTiltPol(ic,inds_list,mp);  
        end

        %--Apply the spectral weights and add together
        Imean = 0; %--ImeanLikeTotallyRight
        for ic=1:Nvals  
            ilam = inds_list(1,ic);
            itt = inds_list(2,ic);
            Imean = Imean + mp.full.lambda_weights_all(ilam)*mp.full.wsTT(itt)/Npol*Iall{ic};  
        end
    else %--Not in parallel
        for si = 1:mp.Nsbp    
            Imean = Imean +  mp.sbp_weights(si)*falco_get_sbp_image(mp,si);
        end
    end

end %--END OF FUNCTION


function Iout = falco_get_single_sim_image_TipTiltPol(ic,inds_list,mp)

ilam = inds_list(1,ic);
itt  = inds_list(2,ic);
ipol = inds_list(3,ic);

%--Get the starlight image
modvar.sbpIndex   = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),1);
modvar.wpsbpIndex = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),2);
mp.full.polaxis = mp.full.pol_conds(ipol);

modvar.whichSource = 'offaxis';
modvar.x_offset = mp.full.xsTT(itt); % used for FALCO full models [lambda0/D]
modvar.y_offset = mp.full.ysTT(itt); % used for FALCO full models [lambda0/D]
mp.full.source_x_offset = mp.full.xsTT(itt); % used for PROPER full models [lambda0/D]
mp.full.source_y_offset = mp.full.ysTT(itt); % used for PROPER full models [lambda0/D]

Estar = model_full(mp, modvar);
Iout = (abs(Estar).^2); %--Apply spectral weighting outside this function

% %--Optionally include the planet PSF
% if(mp.planetFlag)
%     modvar.whichSource = 'exoplanet';
%     Eplanet = model_full(mp,modvar);
%     Iout = Iout + abs(Eplanet).^2; %--Apply spectral weighting outside this function
% end
    
end %--END OF FUNCTION