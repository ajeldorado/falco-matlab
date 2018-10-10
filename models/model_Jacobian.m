% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function jac = model_Jacobian(mp, DM, modvar)
%--Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%
% REVISION HISTORY:
% --------------
% Modified on 2017-11-13 by A.J. Riggs to be compatible with parfor.
% Modified on 2017-11-09 by A.J. Riggs from model_compact.m to
%   model_Jacobian.m
% Modified on 2017-10-17 by A.J. Riggs to have model_compact.m be a wrapper. All the 
%  actual compact models have been moved to sub-routines for clarity.
% Modified on 19 June 2017 by A.J. Riggs to use lower resolution than the
%   full model.
% model_compact.m - 18 August 2016: Modified from hcil_model.m
% hcil_model.m - 18 Feb 2015: Modified from HCIL_model_lab_BB_v3.m
% ---------------
%
% INPUTS:
% -mp = structure of model parameters
% -DM = structure of DM settings
%
% OUTPUTS:
% -jacStruct = control Jacobian for the modes specified

function jacStruct = model_Jacobian(mp)

    %--Calculate the starting DM surfaces beforehand to save time.
    if(any(mp.dm_ind==1)); mp.dm1.compact.surfM = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx,mp.dm1.compact.NdmPad); else; mp.dm1.compact.surfM = zeros(2); end
    if(any(mp.dm_ind==2)); mp.dm2.compact.surfM = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx,mp.dm2.compact.NdmPad); else; mp.dm2.compact.surfM = zeros(2); end
    switch mp.coro
        case{'EHLC'}
            [mp.FPMcube,mp.dm8.surf,mp.dm9.surf] = falco_gen_EHLC_FPM_complex_trans_cube(mp,'compact'); %--Generate the DM surface and FPM outside the full model so that they don't have to be re-made for each wavelength, tip/tilt offset, etc.
            %--Get rid of the DM.dmX.inf_datacube fields in the full model to save RAM.
            mp.dm8 = rmfield(mp.dm8,'inf_datacube'); 
            mp.dm9 = rmfield(mp.dm9,'inf_datacube'); 
        case{'HLC','APHLC','SPHLC'}
            [mp.FPMcube,mp.dm8.surf,mp.dm9.surf] = falco_gen_HLC_FPM_complex_trans_cube(mp,'compact'); %--Generate the DM surface and FPM outside the full model so that they don't have to be re-made for each wavelength, tip/tilt offset, etc.
            %--Get rid of the DM.dmX.inf_datacube fields in the full model to save RAM.
            mp.dm8 = rmfield(mp.dm8,'inf_datacube'); 
            mp.dm9 = rmfield(mp.dm9,'inf_datacube'); 
    end

    %--Initialize the Jacobian cubes for each DM.
    if(any(mp.dm_ind==1)); jacStruct.G1 = zeros(mp.F4.corr.Npix,mp.dm1.Nele,mp.jac.Nmode);  else;  jacStruct.G1 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM1
    if(any(mp.dm_ind==2)); jacStruct.G2 = zeros(mp.F4.corr.Npix,mp.dm2.Nele,mp.jac.Nmode);  else;  jacStruct.G2 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM2
    if(any(mp.dm_ind==3)); jacStruct.G3 = zeros(mp.F4.corr.Npix,mp.dm3.Nele,mp.jac.Nmode);  else;  jacStruct.G3 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM3
    if(any(mp.dm_ind==4)); jacStruct.G4 = zeros(mp.F4.corr.Npix,mp.dm4.Nele,mp.jac.Nmode);  else;  jacStruct.G4 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM4
    if(any(mp.dm_ind==5)); jacStruct.G5 = zeros(mp.F4.corr.Npix,mp.dm5.Nele,mp.jac.Nmode);  else;  jacStruct.G5 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM5
    if(any(mp.dm_ind==6)); jacStruct.G6 = zeros(mp.F4.corr.Npix,mp.dm6.Nele,mp.jac.Nmode);  else;  jacStruct.G6 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM6
    if(any(mp.dm_ind==7)); jacStruct.G7 = zeros(mp.F4.corr.Npix,mp.dm7.Nele,mp.jac.Nmode);  else;  jacStruct.G7 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM7
    if(any(mp.dm_ind==8)); jacStruct.G8 = zeros(mp.F4.corr.Npix,mp.dm8.Nele,mp.jac.Nmode);  else;  jacStruct.G8 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM8
    if(any(mp.dm_ind==9)); jacStruct.G9 = zeros(mp.F4.corr.Npix,mp.dm9.Nele,mp.jac.Nmode);  else;  jacStruct.G9 = zeros(0,0,mp.jac.Nmode);  end % control Jacobian for DM9

    %--Loop over the possible combinations of 1) tip/tilt-offsets, 2) sub-bandpasses, and 3) DM number 
    %   (either with parfor or for)
    fprintf('Computing control Jacobian matrices ... \n'); tic
    vals_list = allcomb(1:mp.jac.Nmode,mp.dm_ind).'; %--dimensions: [2 x length(mp.jac.Nmode)*length(mp.dm_ind) ]
    Nvals = size(vals_list,2);

    
    %--Parallel/distributed computing
    if(mp.flagParfor) 
        parfor ii=1:Nvals
            Jtemp{ii} = model_Jacobian_middle_layer(mp, vals_list, ii);
        end
        
        %--Re-organize the structure
        for ii=1:Nvals
            im = vals_list(1,ii); %--index for tip-tilt-wavelength mode
            whichDM = vals_list(2,ii); %--number of the specified DM

            if(whichDM==1); jacStruct.G1(:,:,im) =  Jtemp{ii};  end
            if(whichDM==2); jacStruct.G2(:,:,im) =  Jtemp{ii};  end
            if(whichDM==8); jacStruct.G8(:,:,im) =  Jtemp{ii};  end
            if(whichDM==9); jacStruct.G9(:,:,im) =  Jtemp{ii};  end
        end
        clear Jtemp

    %--Regular computing (avoid using variable Jtemp to save RAM
    else 
        % fprintf('Jacobian Mode: ');
        for ii=1:Nvals
            im = vals_list(1,ii); %--index for tip-tilt-wavelength mode
            whichDM = vals_list(2,ii); %--number of the specified DM
            fprintf('mode%ddm%d ',im,whichDM)
                        
            if(whichDM==1); jacStruct.G1(:,:,im) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==2); jacStruct.G2(:,:,im) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==8); jacStruct.G8(:,:,im) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==9); jacStruct.G9(:,:,im) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
        end
        fprintf('\n')
    end    

    fprintf('...done.  Time = %.2f\n',toc);

end %--END OF FUNCTION model_Jacobian.m








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function model_Jacobian_middle_layer.m exists at all so that parfor
% can use a linear indexing scheme from 1 to Nmodes.
% I placed the function model_Jacobian_middle_layer.m here as a nested
% function in order to save RAM since the output of the Jacobian structure
% is large and I do not want it copied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% -mp = structure of model parameters
% -vals_list = structure containing combinations of:
%     -tsi = index of the pair of sub-bandpass index and tip/tilt offset index
%     -whichDM = which DM number
%
% OUTPUTS:
% -jacMode = Jacobian for the specified combo of DM, wavelength, and Zernike mode.

function jacMode = model_Jacobian_middle_layer(mp,   vals_list,ii)

    im = vals_list(1,ii); %--index for Zernike-&-subbandpass pair
    whichDM = vals_list(2,ii); %--number of the specified DM

    switch mp.coro 
        case{'EHLC'} %--Extended HLC: DMs, extended FPM with nickel and dielectric modulation, and LS.
            jacMode = model_Jacobian_EHLC(mp,   im, whichDM); 
        
        case{'HLC','APHLC'} %--DMs, optional apodizer, FPM with phase modulation, and LS.
            jacMode = model_Jacobian_HLC(mp,   im, whichDM); 
            
        case{'SPHLC','FHLC'}  %--DMs, optional apodizer, complex/hybrid FPM with outer diaphragm, LS
            jacMode  = model_Jacobian_SPHLC(mp,   im, whichDM); 
            
            
        case{'LC','DMLC','APLC'} %--DMs, optional apodizer, occulting spot FPM, and LS.
            jacMode = model_Jacobian_LC(mp,   im, whichDM); 
            
        case{'SPLC','FLC'} %--DMs, optional apodizer, binary-amplitude FPM with outer diaphragm, LS
            jacMode  = model_Jacobian_SPLC(mp,   im, whichDM); 
            
        case{'vortex','Vortex','VC','AVC'} %--DMs, optional apodizer, vortex FPM, LS
            jacMode  = model_Jacobian_VC(mp,   im, whichDM); 
            
        %case{'SPC','APP','APC'} %--Pupil-plane apodizer is only coronagraphic mask
            %Jac  = model_Jacobian_APC(mp,   tsi, whichDM); 
            
        otherwise
            error('model_Jacobian_middle_layer: CASE NOT RECOGNIZED IN model_Jacobian.m');        
    end    

end %--END OF FUNCTION model_Jacobian_middle_layer.m






