% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% jacStruct : structure containing a separate control Jacobian for each DM

function jacStruct = model_Jacobian(mp)

    %--Calculate the starting DM surfaces beforehand to save time.
    if(any(mp.dm_ind==1)); mp.dm1.compact.surfM = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.NdmPad); else; mp.dm1.compact.surfM = zeros(2); end
    if(any(mp.dm_ind==2)); mp.dm2.compact.surfM = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.NdmPad); else; mp.dm2.compact.surfM = zeros(2); end
    switch upper(mp.coro)
        case{'HLC'}
            if strcmpi(mp.layout, 'fourier')
                [mp.FPMcube, mp.dm8.surf, mp.dm9.surf] = falco_gen_HLC_FPM_complex_trans_cube(mp, 'compact'); %--Generate the DM surface and FPM outside the full model so that they don't have to be re-made for each wavelength, tip/tilt offset, etc.
                %--Get rid of the DM.dmX.inf_datacube fields in the full model to save RAM.
                mp.dm8 = rmfield(mp.dm8, 'inf_datacube'); 
                mp.dm9 = rmfield(mp.dm9, 'inf_datacube'); 
            end
        case{'EHLC'}
            [mp.FPMcube, mp.dm8.surf, mp.dm9.surf] = falco_gen_EHLC_FPM_complex_trans_cube(mp, 'compact'); %--Generate the DM surface and FPM outside the full model so that they don't have to be re-made for each wavelength, tip/tilt offset, etc.
            %--Get rid of the DM.dmX.inf_datacube fields in the full model to save RAM.
            mp.dm8 = rmfield(mp.dm8, 'inf_datacube'); 
            mp.dm9 = rmfield(mp.dm9, 'inf_datacube'); 
    end

    %--Initialize the Jacobian cubes for each DM.
    if(mp.flagFiber)
        if(mp.flagLenslet)
            if(any(mp.dm_ind==1)); jacStruct.G1 = zeros(mp.Fend.Nlens, mp.dm1.Nele, mp.jac.Nmode);  else;  jacStruct.G1 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM1
            if(any(mp.dm_ind==2)); jacStruct.G2 = zeros(mp.Fend.Nlens, mp.dm2.Nele, mp.jac.Nmode);  else;  jacStruct.G2 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM2
        else
        	if(any(mp.dm_ind==1)); jacStruct.G1 = zeros(mp.Fend.corr.Npix, mp.dm1.Nele, mp.jac.Nmode);  else;  jacStruct.G1 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM1
            if(any(mp.dm_ind==2)); jacStruct.G2 = zeros(mp.Fend.corr.Npix, mp.dm2.Nele, mp.jac.Nmode);  else;  jacStruct.G2 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM2
        end
        jacStruct.G3 = zeros(0, 0, mp.jac.Nmode);
        jacStruct.G4 = zeros(0, 0, mp.jac.Nmode);
        jacStruct.G5 = zeros(0, 0, mp.jac.Nmode);
        jacStruct.G6 = zeros(0, 0, mp.jac.Nmode);
        jacStruct.G7 = zeros(0, 0, mp.jac.Nmode);
        jacStruct.G8 = zeros(0, 0, mp.jac.Nmode);
        jacStruct.G9 = zeros(0, 0, mp.jac.Nmode);
    else
        if(any(mp.dm_ind==1)); jacStruct.G1 = zeros(mp.Fend.corr.Npix, mp.dm1.Nele, mp.jac.Nmode);  else;  jacStruct.G1 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM1
        if(any(mp.dm_ind==2)); jacStruct.G2 = zeros(mp.Fend.corr.Npix, mp.dm2.Nele, mp.jac.Nmode);  else;  jacStruct.G2 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM2
        if(any(mp.dm_ind==3)); jacStruct.G3 = zeros(mp.Fend.corr.Npix, mp.dm3.Nele, mp.jac.Nmode);  else;  jacStruct.G3 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM3
        if(any(mp.dm_ind==4)); jacStruct.G4 = zeros(mp.Fend.corr.Npix, mp.dm4.Nele, mp.jac.Nmode);  else;  jacStruct.G4 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM4
        if(any(mp.dm_ind==5)); jacStruct.G5 = zeros(mp.Fend.corr.Npix, mp.dm5.Nele, mp.jac.Nmode);  else;  jacStruct.G5 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM5
        if(any(mp.dm_ind==6)); jacStruct.G6 = zeros(mp.Fend.corr.Npix, mp.dm6.Nele, mp.jac.Nmode);  else;  jacStruct.G6 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM6
        if(any(mp.dm_ind==7)); jacStruct.G7 = zeros(mp.Fend.corr.Npix, mp.dm7.Nele, mp.jac.Nmode);  else;  jacStruct.G7 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM7
        if(any(mp.dm_ind==8)); jacStruct.G8 = zeros(mp.Fend.corr.Npix, mp.dm8.Nele, mp.jac.Nmode);  else;  jacStruct.G8 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM8
        if(any(mp.dm_ind==9)); jacStruct.G9 = zeros(mp.Fend.corr.Npix, mp.dm9.Nele, mp.jac.Nmode);  else;  jacStruct.G9 = zeros(0, 0, mp.jac.Nmode);  end % control Jacobian for DM9
    end
        
    %--Loop over the possible combinations of 1) Zernike modes, 2) sub-bandpasses, 3) stars, & 4) DM number 
    %   (either with parfor or for)
    fprintf('Computing control Jacobian matrices ... \n'); tic
    vals_list = allcomb(1:mp.jac.Nmode, mp.dm_ind).'; %--dimensions: [2 x length(mp.jac.Nmode)*length(mp.dm_ind) ]
    Nvals = size(vals_list, 2);

    %--Parallel/distributed computing
    if(mp.flagParfor)
        if isfield(mp, 'tb')
            mp = rmfield(mp, 'tb'); % Remove testbed object 'tb' if it exists before calling parfor.
        end
        parfor ii=1:Nvals
            Jtemp{ii} = model_Jacobian_middle_layer(mp, vals_list, ii)
        end        
        
        %--Re-organize the structure
        for ii = 1:Nvals
            iMode = vals_list(1, ii); % index for Zernike-wavelength-star mode
            whichDM = vals_list(2, ii); % number of the specified DM

            if(whichDM==1); jacStruct.G1(:, :, iMode) =  Jtemp{ii};  end
            if(whichDM==2); jacStruct.G2(:, :, iMode) =  Jtemp{ii};  end
            if(whichDM==3); jacStruct.G3(:, :, iMode) =  Jtemp{ii};  end
            if(whichDM==4); jacStruct.G4(:, :, iMode) =  Jtemp{ii};  end
            if(whichDM==5); jacStruct.G5(:, :, iMode) =  Jtemp{ii};  end
            if(whichDM==6); jacStruct.G6(:, :, iMode) =  Jtemp{ii};  end
            if(whichDM==7); jacStruct.G7(:, :, iMode) =  Jtemp{ii};  end
            if(whichDM==8); jacStruct.G8(:, :, iMode) =  Jtemp{ii};  end
            if(whichDM==9); jacStruct.G9(:, :, iMode) =  Jtemp{ii};  end
        end
        clear Jtemp

    %--Serial calculation
    else 
        for ii = 1:Nvals
            iMode = vals_list(1, ii); %--index for Zernike-wavelength-star mode
            whichDM = vals_list(2, ii); %--number of the specified DM
            fprintf('mode%ddm%d ', iMode, whichDM)
                        
            if(whichDM==1); jacStruct.G1(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==2); jacStruct.G2(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==3); jacStruct.G3(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==4); jacStruct.G4(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==5); jacStruct.G5(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==6); jacStruct.G6(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==7); jacStruct.G7(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==8); jacStruct.G8(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
            if(whichDM==9); jacStruct.G9(:, :, iMode) =  model_Jacobian_middle_layer(mp, vals_list, ii);  end
        end
        fprintf('\n')
    end    

    fprintf('...done.  Time = %.2f\n', toc);

    
    %% TIED ACTUATORS
    
    %--Handle tied actuators by adding the 2nd actuator's Jacobian column to
    %the first actuator's column, and then zeroing out the 2nd actuator's column.
    if(any(mp.dm_ind==1))
        mp.dm1 = falco_enforce_dm_constraints(mp.dm1); %-Update the sets of tied actuators
        for ti=1:size(mp.dm1.tied, 1)
            Index1all = mp.dm1.tied(ti, 1); %--Index of first tied actuator in whole actuator set. 
            Index2all = mp.dm1.tied(ti, 2); %--Index of second tied actuator in whole actuator set. 
            Index1subset = find(mp.dm1.act_ele==Index1all); %--Index of first tied actuator in subset of used actuators. 
            Index2subset = find(mp.dm1.act_ele==Index2all); %--Index of second tied actuator in subset of used actuators. 
            jacStruct.G1(:, Index1subset, :) = jacStruct.G1(:, Index1subset, :) + jacStruct.G1(:, Index2subset, :); % adding the 2nd actuators Jacobian column to the first actuator's column
            jacStruct.G1(:, Index2subset, :) = 0*jacStruct.G1(:, Index2subset, :); % zero out the 2nd actuator's column.
        end
    end
    if(any(mp.dm_ind==2))
        mp.dm2 = falco_enforce_dm_constraints(mp.dm2); %-Update the sets of tied actuators
        for ti=1:size(mp.dm2.tied, 1)
            Index1all = mp.dm2.tied(ti, 1); %--Index of first tied actuator in whole actuator set. 
            Index2all = mp.dm2.tied(ti, 2); %--Index of second tied actuator in whole actuator set. 
            Index1subset = find(mp.dm2.act_ele==Index1all); %--Index of first tied actuator in subset of used actuators. 
            Index2subset = find(mp.dm2.act_ele==Index2all); %--Index of second tied actuator in subset of used actuators. 
            jacStruct.G2(:, Index1subset, :) = jacStruct.G2(:, Index1subset, :) + jacStruct.G2(:, Index2subset, :); % adding the 2nd actuators Jacobian column to the first actuator's column
            jacStruct.G2(:, Index2subset, :) = 0*jacStruct.G2(:, Index2subset, :); % zero out the 2nd actuator's column.
        end
    end
    if(any(mp.dm_ind==8))
        for ti=1:size(mp.dm8.tied, 1)
            Index1all = mp.dm8.tied(ti, 1); %--Index of first tied actuator in whole actuator set. 
            Index2all = mp.dm8.tied(ti, 2); %--Index of second tied actuator in whole actuator set. 
            Index1subset = find(mp.dm8.act_ele==Index1all); %--Index of first tied actuator in subset of used actuators. 
            Index2subset = find(mp.dm8.act_ele==Index2all); %--Index of second tied actuator in subset of used actuators. 
            jacStruct.G8(:, Index1subset, :) = jacStruct.G8(:, Index1subset, :) + jacStruct.G8(:, Index2subset, :); % adding the 2nd actuators Jacobian column to the first actuator's column
            jacStruct.G8(:, Index2subset, :) = 0*jacStruct.G8(:, Index2subset, :); % zero out the 2nd actuator's column.
        end
    end
    if(any(mp.dm_ind==9))
        for ti=1:size(mp.dm9.tied, 1)
            Index1all = mp.dm9.tied(ti, 1); %--Index of first tied actuator in whole actuator set. 
            Index2all = mp.dm9.tied(ti, 2); %--Index of second tied actuator in whole actuator set. 
            Index1subset = find(mp.dm9.act_ele==Index1all); %--Index of first tied actuator in subset of used actuators. 
            Index2subset = find(mp.dm9.act_ele==Index2all); %--Index of second tied actuator in subset of used actuators. 
            jacStruct.G9(:, Index1subset, :) = jacStruct.G9(:, Index1subset, :) + jacStruct.G9(:, Index2subset, :); % adding the 2nd actuators Jacobian column to the first actuator's column
            jacStruct.G9(:, Index2subset, :) = 0*jacStruct.G9(:, Index2subset, :); % zero out the 2nd actuator's column.
        end
    end
       
end %--END OF FUNCTION model_Jacobian.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select which optical layout's Jacobian model to use and get the output E-field
% INPUTS
% ------
% mp : structure of model parameters
% vals_list : structure containing combinations of:
%     -tsi : index of the pair of sub-bandpass index and tip/tilt offset index
%     -whichDM : which DM number
%
% OUTPUTS
% -------
% -jacMode = Jacobian for the specified combo of DM, wavelength, and Zernike mode.
%
function jacMode = model_Jacobian_middle_layer(mp, vals_list, jj)

    iModeCopy = vals_list(1, jj); %--index for Zernike-&-subbandpass pair
    whichDMCopy = vals_list(2, jj); %--number of the specified DM

    switch upper(mp.coro)
        
        case{'LC', 'APLC', 'SPLC', 'FLC', 'HLC'}
            jacMode = model_Jacobian_lyot(mp, iModeCopy, whichDMCopy); 
            
        case{'VORTEX', 'VC'}
            jacMode = model_Jacobian_VC(mp, iModeCopy, whichDMCopy);
            
        case{'EHLC'} %--Extended HLC: DMs, extended FPM with nickel and dielectric modulation, and LS.
            jacMode = model_Jacobian_EHLC(mp, iModeCopy, whichDMCopy);
            
        otherwise
            error('No Jacobian function for the value of mp.coro');        
    end

end %--END OF FUNCTION model_Jacobian_middle_layer.m    
