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

function [Emat,Isum2D] = falco_est_perfect_Efield_compact(mp);
    
    if(isfield(mp,'lowfs')==false)
        mp.lowfs = false; %--Set LOWFS flag to false if it isn't included
    end

    IfocusCube = zeros(mp.Fend.Neta, mp.Fend.Nxi, mp.Nsbp);
    Emat = zeros(mp.Fend.corr.Npix, mp.Nsbp);
    
    for si=1:mp.Nsbp
        modvar.sbpIndex = si;
        modvar.zernIndex = 1;
        modvar.wpsbpIndex = mp.wi_ref;
        modvar.whichSource = 'star';

        E2D = model_compact(mp, modvar);

        if mp.lowfs
            Icube(:,:,si) = abs(E2D).^2;
        else
            Emat(:,si) = E2D(mp.Fend.corr.inds); % Actual field in estimation area
            IfocusCube(:,:,si) = (abs(E2D).^2)*mp.jac.weightMat(si,1);
        end

    end
    
    if mp.lowfs
        Isum2D = Icube; %--Still want the cube for the LOWFS
    else
        Isum2D = sum(IfocusCube,3);
    end

end %--END OF FUNCTION