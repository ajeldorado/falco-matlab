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
        if(mp.lowfs)
            error('falco_est_perfect_Efield_with_Zernikes.m: Do not call this function to retrieve the LOWFS E-field. ')
        end
    end
   
    if(mp.flagFiber)
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
        if(mp.flagFiber)
            [I, J] = ind2sub(size(mp.F5.RHOS), find(~mp.F5.RHOS));
            Emat(:,im) = EfibCompact(I,J,:);
        else
            Emat(:,im) = E2D(mp.Fend.corr.inds); % Actual field in estimation area
        end
    end
    
end %--END OF FUNCTION