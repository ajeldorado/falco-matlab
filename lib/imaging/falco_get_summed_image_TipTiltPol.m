% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Get a broadband image over the entire bandpass by summing over subbands,
% tip/tilt settings, and polarization states.
%
% INPUTS
% ------
% mp : structure of all model parameters
%
% OUTPUTS
% -------
% Itotal : image in units of normalized intensity

function Itotal = falco_get_summed_image_TipTiltPol(mp)

    %--Compute the DM surfaces outside the full model to save some time
    if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
    if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
    if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end

    Itotal = 0; % Initialize image
    
    %--Number of polarization states used
    mp.full.dummy = 1; %--Initialize if this doesn't exist
    if(isfield(mp.full,'pol_conds'))  
        Npol = length(mp.full.pol_conds);  
    else
        Npol = 1;
    end
    
    %--Generate the tip/tilt offsets and their normalized weights
    [mp.full.xsTT, mp.full.ysTT, mp.full.wsTT] = falco_gen_RMS_TipTilt(mp.full.TTrms,mp.full.Dstar,mp.full.Dtel,mp.lambda0,'Nacross',mp.full.TipTiltNacross);
    Ntt = length(mp.full.wsTT);
    fprintf('\n%d tip-tilt offset points used.\n',Ntt);
    %--Iterate over all combinations of sub-bandpass, wavelength, tip/tilt offset, and polarization state.
    
    if(mp.flagSim)
        %--Loop over all wavelengths, tip/tilt offsets, and polarizations        
        inds_list = allcomb(1:mp.full.NlamUnique,1:Ntt,1:Npol).'; %--dimensions: [3 x mp.full.NlamUnique*Ntt*Npol ]
        %inds_list = allcomb(1:mp.Nsbp,1:mp.Nwpsbp,1:Ntt,1:Npol).'; %--dimensions: [4 x mp.Nsbp*Nwpsbp*Ntt*Npol ]
        Nvals = size(inds_list,2);
        
        if(mp.flagParfor) %--Save a lot of time by running PROPER full model in parallel
            parfor ic=1:Nvals
                Iall{ic} = falco_get_single_sim_image_TipTiltPol(ic,inds_list,mp);  
            end
        else
            for ic=Nvals:-1:1
                Iall{ic} = falco_get_single_sim_image_TipTiltPol(ic,inds_list,mp);  
            end
        end

        %--Apply the spectral weights and add together
        Itotal = 0; %--ItotalLikeTotallyRight
        for ic=1:Nvals  
            ilam = inds_list(1,ic);
            itt = inds_list(2,ic);
            Itotal = Itotal + mp.full.lambda_weights_all(ilam)*mp.full.wsTT(itt)/Npol*Iall{ic};  
        end
    end

end %--END OF FUNCTION


function Iout = falco_get_single_sim_image_TipTiltPol(ic,inds_list,mp)

ilam = inds_list(1,ic);
itt  = inds_list(2,ic);
ipol = inds_list(3,ic);
modvar.starIndex = 1;

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
    
end %--END OF FUNCTION