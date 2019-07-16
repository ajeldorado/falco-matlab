% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to return the perfect-knowledge E-field for the compact model at
% each wavelength and 1st-order Zernike mode.
%
% REVISION HISTORY:
%--Modified on 2018-09-24 by A.J. Riggs to include the E-fields from the
%    1st-order Zernike modes.
%--Modified on 2018-08-27 by A.J. Riggs to include the LOWFS parts from
%    Erkin's code.
%--Created on 2018-01-24 by A.J. Riggs.

function Emat = falco_est_perfect_Efield_with_Zernikes(mp)

    if(isfield(mp,'lowfs'))
        if(mp.lowfs);  error('falco_est_perfect_Efield_with_Zernikes.m: Do not call this function to retrieve the LOWFS E-field. ');  end
    end

    %--Polarization states as defined in the guide by John Krist for the
    % WFIRST CGI Phase B model.
    if(isfield(mp.full,'pol_conds'))
        if( isempty(setdiff([-2,-1,1,2],mp.full.pol_conds)) ) %--X and Y out
            mp.full.polaxis = 10; %--Use the average polarization state for the perfect estimation.
        elseif( isempty(setdiff([-1,1],mp.full.pol_conds)) ) %--X out
            mp.full.polaxis = 5; %--Use the average polarization state for the perfect estimation.
        elseif( isempty(setdiff([-2,2],mp.full.pol_conds)) ) %--Y out
            mp.full.polaxis = 6; %--Use the average polarization state for the perfect estimation.    
        end
    end

    if(mp.flagFiber)
        if(mp.flagLenslet)
            Emat = zeros(mp.Fend.Nlens, mp.jac.Nmode);
        else
            Emat = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
        end
        
        for im=1:mp.jac.Nmode
            modvar.sbpIndex = mp.jac.sbp_inds(im);
            modvar.zernIndex = mp.jac.zern_inds(im);
            modvar.wpsbpIndex = mp.wi_ref;
            modvar.whichSource = 'star';

            [E2D, EfibCompact] = model_full(mp, modvar);
            
            if(mp.flagLenslet)
                [I, J] = ind2sub(size(mp.F5.RHOS), find(~mp.F5.RHOS));
                Emat(:,im) = EfibCompact(I,J,:);
            else
                Emat(:,im) = EfibCompact(mp.Fend.corr.inds);
            end
        end

    else %--Regular, non-fiber simulations
        Emat = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
        
        if(mp.flagParfor) %--Perform all cases in parallel
            %--Loop over all modes and wavelengths       
            ind_list = allcomb(1:mp.jac.Nmode,1:mp.Nwpsbp).';
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
            
            counter = 0;
            for im=1:mp.jac.Nmode
                EsbpMean = 0;
                for wi=1:mp.Nwpsbp
                    counter = counter + 1;
                    EsbpMean = EsbpMean + EmatAll(:,counter)*mp.full.lambda_weights(wi);
                end
                Emat(:,im) = EsbpMean;
            end
            
        else %--Not done in parallel
            for im=1:mp.jac.Nmode
                modvar.sbpIndex = mp.jac.sbp_inds(im);
                modvar.zernIndex = mp.jac.zern_inds(im);
                modvar.whichSource = 'star';
                
                %--Take the mean over the wavelengths within the sub-bandpass
                EmatSbp = zeros(mp.Fend.corr.Npix, mp.Nwpsbp);
                for wi=1:mp.Nwpsbp
                    modvar.wpsbpIndex = wi;
                    E2D = model_full(mp, modvar);
                    EmatSbp(:,wi) = mp.full.lambda_weights(wi)*E2D(mp.Fend.corr.inds); %--Actual field in estimation area. Apply spectral weight within the sub-bandpass
                end
                Emat(:,im) = sum(EmatSbp,2);
            end
        end        
    end
    
end %--END OF FUNCTION

%--Extra function needed to use parfor (because parfor can have only a
%  single changing input argument).
function Evec = falco_est_perfect_Efield_with_Zernikes_parfor(ni,ind_list,mp)

    im = ind_list(1,ni); %--Index of the Jacobian mode
    wi = ind_list(2,ni); %--Index of the wavelength in the sub-bandpass
    
    modvar.sbpIndex = mp.jac.sbp_inds(im);
    modvar.zernIndex = mp.jac.zern_inds(im);
    modvar.wpsbpIndex = wi;
    modvar.whichSource = 'star';
    
    E2D = model_full(mp, modvar);
    Evec = E2D(mp.Fend.corr.inds); % Actual field in estimation area. Don't apply spectral weight here.

end

