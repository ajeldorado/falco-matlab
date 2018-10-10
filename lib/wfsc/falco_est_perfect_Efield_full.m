% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to return the perfect-knowledge E-field and summed intensity for
% the full model.
%
% REVISION HISTORY: 
% -Created on 2018-01-24 by A.J. Riggs.

function [Emat,Imat] = falco_est_perfect_Efield_full(mp,DM)
    
%     Icube = zeros(mp.F4.Neta, mp.F4.Nxi, mp.jac.Nmode);
    Imat = zeros(mp.F4.Neta, mp.F4.Nxi);
    Emat = zeros(mp.F4.corr.Npix, mp.jac.Nmode);
    
    for im=1:mp.jac.Nmode
        % Get full-knowledge, aberrated, noiseless starlight E-field and image
        modvar.flagCalcJac = 0; 
        modvar.sbpIndex = mp.jac.sbp_inds(im); %mp.Wttlam_si(tsi);
        modvar.zernIndex = mp.jac.zern_inds(im);
        %modvar.ttIndex = mp.jac.zern_inds(im);
        modvar.wpsbpIndex = mp.wi_ref;
        %modvar.lambda = mp.sbp_centers(modvar.sbpIndex);
        modvar.whichSource = 'star';

        E2D = model_full(mp, DM, modvar);
        Emat(:,im) = E2D(mp.F4.corr.maskBool);   % Actual field in estimation area
        Imat = Imat + (abs(E2D).^2)*mp.jac.weights(im);
        %Icube(:,:,im) = (abs(E2D).^2)*mp.jac.weights(im);
    end
    %Imat = sum(Icube,3);


end %--END OF FUNCTION
