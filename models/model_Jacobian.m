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

function jacStruct = model_Jacobian(mp, DM)

    %--Calculate the starting DM surfaces beforehand to save time.
    %--Compute the DM surfaces outside the full model to save lots of time
    if(any(DM.dm_ind==1)); DM.dm1.compact.surfM = falco_gen_dm_surf(DM.dm1, DM.dm1.compact.dx,DM.dm1.compact.NdmPad); else DM.dm1.compact.surfM = zeros(2); end;
    if(any(DM.dm_ind==2)); DM.dm2.compact.surfM = falco_gen_dm_surf(DM.dm2, DM.dm2.compact.dx,DM.dm2.compact.NdmPad); else DM.dm2.compact.surfM = zeros(2); end;


    %--Get rid of the DM.dmX.inf_datacube fields in the full model because they are HUGE and will
    % take up way too much RAM if passed to all threads. Just pass in the one
    % influence function that will be used.
    if(any(DM.dm_ind==1)); G1 = zeros(length(mp.F4.compact.corr.inds),DM.dm1.NactTotal,mp.Nttlam); end % control Jacobian for DM1
    if(any(DM.dm_ind==2)); G2 = zeros(length(mp.F4.compact.corr.inds),DM.dm2.NactTotal,mp.Nttlam); end % control Jacobian for DM2


    %--Loop over the possible combinations of 1) tip/tilt-offsets, 2) sub-bandpasses, and 3) DM number 
    %   (either with parfor or for)
    fprintf('Computing control Jacobian matrices ... \n'); tic
    vals_list = allcomb(1:mp.Nttlam,DM.dm_ind); %--dimensions: [2 x length(mp.Nttlam)*length(DM.dm_ind) ]
    Nvals = size(vals_list,2);

    %--Parallel/distributed computing
    if(mp.flagParfor) 
        parfor ii=1:Nvals
            Jtemp{ii} = model_Jacobian_middle_layer(mp, DM, vals_list, ii);
        end

    %--Regular computing    
    else 
        % fprintf('Jacobian Mode: ');
        for ii=1:Nvals
            tsi = vals_list(1,ii); %--index for tip-tilt-wavelength mode
            whichDM = vals_list(2,ii); %--number of the specified DM
            fprintf('Mode%dDM%d ',tsi,whichDM)
            Jtemp{ii} = model_Jacobian_middle_layer(mp, DM, vals_list, ii);
        end
        fprintf('\n')
    end    


    %--Re-organize the structure
    for ii=1:Nvals
        tsi = vals_list(1,ii); %--index for tip-tilt-wavelength mode
        whichDM = vals_list(2,ii); %--number of the specified DM

        if(whichDM==1); G1(:,:,tsi) =  Jtemp{ii};  end
        if(whichDM==2); G2(:,:,tsi) =  Jtemp{ii};  end
    end
    clear Jtemp

    if(any(DM.dm_ind==1)); jacStruct.G1 = G1; end
    if(any(DM.dm_ind==2)); jacStruct.G2 = G2; end

    fprintf('...done.  Time = %.2f\n',toc);

end % End of function
