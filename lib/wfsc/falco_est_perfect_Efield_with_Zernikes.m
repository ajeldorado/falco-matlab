% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to return the perfect-knowledge E-field from the full model for
% the all the modes (combinations of subband, Zernike, and star)
%
% INPUTS
% mp : structure of model parameters
%
% OUTPUTS
% Emat : 2-D complex-valued array of E-field values. Dimensions are number
% of pixels in the dark hole by the number of Jacobian modes.


function Emat = falco_est_perfect_Efield_with_Zernikes(mp)

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
            Emat = zeros(mp.Fend.Nlens, mp.jac.Nmode);
        else
            Emat = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
        end
        
        for iMode=1:mp.jac.Nmode
            modvar.sbpIndex = mp.jac.sbp_inds(iMode);
            modvar.zernIndex = mp.jac.zern_inds(iMode);
            modvar.wpsbpIndex = mp.wi_ref;
            modvar.whichSource = 'star';

            [E2D, EfiberCompact] = model_full(mp, modvar);
            
            if(mp.flagLenslet)
                [I, J] = ind2sub(size(mp.F5.RHOS), find(~mp.F5.RHOS));
                Emat(:, iMode) = EfiberCompact(I, J, :);
            else
                Emat(:, iMode) = EfiberCompact(mp.Fend.corr.inds);
            end
        end

    else %--Regular (non-fiber) simulations
        Emat = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
        
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
                EsbpMean = 0;
                for wi=1:mp.Nwpsbp
                    counter = counter + 1;
                    EsbpMean = EsbpMean + EmatAll(:,counter)*mp.full.lambda_weights(wi);
                end
                Emat(:, iMode) = EsbpMean;
            end
            
        else %--Not done in parallel
            for iMode = 1:mp.jac.Nmode
                modvar.sbpIndex = mp.jac.sbp_inds(iMode);
                modvar.zernIndex = mp.jac.zern_inds(iMode);
                modvar.starIndex = mp.jac.star_inds(iMode);
                modvar.whichSource = 'star';
                
                %--Take the mean over the wavelengths within the sub-bandpass
                EmatSbp = zeros(mp.Fend.corr.Npix, mp.Nwpsbp);
                for wi = 1:mp.Nwpsbp
                    modvar.wpsbpIndex = wi;
                    E2D = model_full(mp, modvar);
                    EmatSbp(:, wi) = mp.full.lambda_weights(wi)*E2D(mp.Fend.corr.inds); 
                end
                Emat(:, iMode) = sum(EmatSbp, 2);
            end
        end        
    end
    
end %--END OF FUNCTION

%--Extra function needed to use parfor (because parfor can have only a
%  single changing input argument).
function Evec = falco_est_perfect_Efield_with_Zernikes_parfor(ni, ind_list, mp)

    iMode = ind_list(1,ni); %--Index of the Jacobian mode
    wi = ind_list(2,ni); %--Index of the wavelength in the sub-bandpass
    
    modvar.sbpIndex = mp.jac.sbp_inds(iMode);
    modvar.zernIndex = mp.jac.zern_inds(iMode);
    modvar.starIndex = mp.jac.star_inds(iMode);
    modvar.wpsbpIndex = wi;
    modvar.whichSource = 'star';
    
    E2D = model_full(mp, modvar);
    Evec = E2D(mp.Fend.corr.inds);
end

