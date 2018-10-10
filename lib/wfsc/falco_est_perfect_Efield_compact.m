% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to return the perfect-knowledge E-field and summed intensity for
% the compact model.
%
% REVISION HISTORY:
%--Modified on 2018-08-27 by A.J. Riggs to include the LOWFS parts from
%    Erkin's code.
%--Created on 2018-01-24 by A.J. Riggs.

function [Emat,Isum2D] = falco_est_perfect_Efield_compact(mp,DM)
    
    if(isfield(mp,'lowfs')==false)
        mp.lowfs = false; %--Set LOWFS flag to false if it isn't included
    end

    
    IfocusCube = zeros(mp.F4.Neta, mp.F4.Nxi, mp.Nsbp);
    Emat = zeros(mp.F4.corr.Npix, mp.Nsbp);
    
    for si=1:mp.Nsbp
        modvar.flagCalcJac = 0; 
        modvar.sbpIndex = si; %mp.jac.sbp_inds(im); %mp.Wttlam_si(im);
        modvar.zernIndex = 1;%mp.jac.zern_inds(im);
        %modvar.ttIndex = mp.Wttlam_ti(im);
        modvar.wpsbpIndex = mp.wi_ref;
        modvar.whichSource = 'star';

        E2D = model_compact(mp, DM, modvar);

        if mp.lowfs
            Icube(:,:,si) = abs(E2D).^2;
        else
            Emat(:,si) = E2D(mp.F4.corr.inds);   % Actual field in estimation area
            IfocusCube(:,:,si) = (abs(E2D).^2)*mp.jac.weightMat(si,1);%mp.jac.weights(im);
        end

    end
    
    
%     IfocusCube = zeros(mp.F4.Neta, mp.F4.Nxi, mp.jac.Nmode);
%     Emat = zeros(mp.F4.corr.Npix, mp.jac.Nmode);
%     
%     for im=1:mp.jac.Nmode
%         modvar.flagCalcJac = 0; 
%         modvar.sbpIndex = mp.jac.sbp_inds(im); %mp.Wttlam_si(im);
%         modvar.zernIndex = mp.jac.zern_inds(im);
%         %modvar.ttIndex = mp.Wttlam_ti(im);
%         modvar.wpsbpIndex = mp.wi_ref;
%         modvar.whichSource = 'star';
% 
%         E2D = model_compact(mp, DM, modvar);
% 
%         if mp.lowfs
%             Icube(:,:,im) = abs(E2D).^2;
%         else
%             Emat(:,im) = E2D(mp.F4.corr.inds);   % Actual field in estimation area
%             IfocusCube(:,:,im) = (abs(E2D).^2)*mp.jac.weights(im);
%         end
% 
%     end
    
    
    if mp.lowfs
        Isum2D = Icube; %--Still want the cube for the LOWFS
    else
        Isum2D = sum(IfocusCube,3);
    end



end %--END OF FUNCTION
