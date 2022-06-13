% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Return the perfect-knowledge E-field from the full model for
% the all the modes (combinations of subband, Zernike, and star)
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% ev : structure of estimation variables

function ev = falco_est_perfect_Efield_with_Zernikes(mp)

    ev.imageArray = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1, mp.Nsbp);

    %--Polarization states as defined in the guide by John Krist for the
    % WFIRST CGI Phase B model.
    if isfield(mp.full, 'pol_conds')
        if isempty(setdiff([-2,-1,1,2], mp.full.pol_conds)) %--X and Y out
            mp.full.polaxis = 10; %--Use the average polarization state for the perfect estimation.
        elseif isempty(setdiff([-1,1], mp.full.pol_conds)) %--X out
            mp.full.polaxis = 5; %--Use the average polarization state for the perfect estimation.
        elseif isempty(setdiff([-2,2], mp.full.pol_conds)) %--Y out
            mp.full.polaxis = 6; %--Use the average polarization state for the perfect estimation.    
        end
    end

    if(mp.flagFiber)
        if(mp.flagLenslet)
            Eest = zeros(mp.Fend.Nlens, mp.jac.Nmode);
        else
            Eest = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
        end
        
        for iMode=1:mp.jac.Nmode
            modvar = ModelVariables;
            modvar.sbpIndex = mp.jac.sbp_inds(iMode);
            modvar.zernIndex = mp.jac.zern_inds(iMode);
            modvar.wpsbpIndex = mp.wi_ref;
            modvar.whichSource = 'star';

            [E2D, EfiberCompact] = model_full(mp, modvar);
            
            if(mp.flagLenslet)
                [I, J] = ind2sub(size(mp.F5.RHOS), find(~mp.F5.RHOS));
                Eest(:, iMode) = EfiberCompact(I, J, :);
            else
                Eest(:, iMode) = EfiberCompact(mp.Fend.corr.inds);
            end
        end

    else %--Regular (non-fiber) simulations
        Eest = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
        
        if(mp.flagParfor) %--Perform all cases in parallel
            %--Loop over all modes and wavelengths       
            ind_list = allcomb(1:mp.jac.Nmode, 1:mp.Nwpsbp).';
            Nval = size(ind_list,2);
            
            parfor ni=1:Nval
                Evecs{ni} = falco_est_perfect_Efield_with_Zernikes_parfor(ni,ind_list,mp)
            end
            
            %--Re-order for easier indexing
            EmatAll = zeros(mp.Fend.corr.Npix, Nval);
            for ni=1:Nval
                EmatAll(:,ni) = Evecs{ni};
            end
            clear Evecs
            
            % Average the E-field over the wavelengths in a subband
            counter = 0;
            for iMode = 1:mp.jac.Nmode
                EsubbandMean = 0;
                for wi=1:mp.Nwpsbp
                    counter = counter + 1;
                    EsubbandMean = EsubbandMean + mp.full.lambda_weights(wi)*EmatAll(:, counter);
                end
                Eest(:, iMode) = EsubbandMean;
                
                iSubband = mp.jac.sbp_inds(iMode);
                iZernike = mp.jac.zern_inds(iMode);
                isPiston = (iZernike == 1);
                if isPiston
                    imageTemp = zeros(mp.Fend.Neta, mp.Fend.Nxi);
                    imageTemp(mp.Fend.corr.maskBool) = abs(EsubbandMean).^2;
                    ev.imageArray(:, :, 1, iSubband) = imageTemp + ev.imageArray(:, :, 1, iSubband);
                end
            end
            
        else %--Not done in parallel
            for iMode = 1:mp.jac.Nmode
                modvar = ModelVariables;
                modvar.sbpIndex = mp.jac.sbp_inds(iMode);
                modvar.zernIndex = mp.jac.zern_inds(iMode);
                modvar.starIndex = mp.jac.star_inds(iMode);
                modvar.whichSource = 'star';
                isPiston = (modvar.zernIndex == 1);
                
                %--Take the mean over the wavelengths within the sub-bandpass
                EmatSubband = zeros(mp.Fend.Neta, mp.Fend.Nxi);
                for wi = 1:mp.Nwpsbp
                    modvar.wpsbpIndex = wi;
                    E2D = model_full(mp, modvar);
                    EmatSubband = EmatSubband + mp.full.lambda_weights(wi)*E2D; 
                end
                if isPiston
                    ev.imageArray(:, :, 1, modvar.sbpIndex) = abs(EmatSubband).^2;
                end
                Eest(:, iMode) = EmatSubband(mp.Fend.corr.maskBool);
                
            end
        end        
    end
    
    ev.Eest = Eest;
    ev.IincoEst = zeros(size(ev.Eest));
    
end %--END OF FUNCTION

%--Extra function needed to use parfor (because parfor can have only a
%  single changing input argument).
function Evec = falco_est_perfect_Efield_with_Zernikes_parfor(ni, ind_list, mp)

    iMode = ind_list(1,ni); %--Index of the Jacobian mode
    wi = ind_list(2,ni); %--Index of the wavelength in the sub-bandpass

    modvar = ModelVariables;
    modvar.sbpIndex = mp.jac.sbp_inds(iMode);
    modvar.zernIndex = mp.jac.zern_inds(iMode);
    modvar.starIndex = mp.jac.star_inds(iMode);
    modvar.wpsbpIndex = wi;
    modvar.whichSource = 'star';
    
    E2D = model_full(mp, modvar);
    Evec = E2D(mp.Fend.corr.maskBool);
end
