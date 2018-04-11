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
% -Created on 2018-01-24 by A.J. Riggs.

function [Emat,Isum2D] = falco_est_perfect_Efield_compact(mp,DM)
    
    Icube = zeros(mp.F4.compact.Neta, mp.F4.compact.Nxi, mp.Nttlam);
    Emat = zeros(length(mp.F4.compact.corr.inds), mp.Nttlam);
    
    for tsi=1:mp.Nttlam
        modvar.flagCalcJac = 0; 
        modvar.sbpIndex = mp.Wttlam_si(tsi);
        modvar.ttIndex = mp.Wttlam_ti(tsi);
        modvar.wpsbpIndex = mp.wi_ref;
        modvar.whichSource = 'star';

        E2D = model_compact(mp, DM, modvar);
        Emat(:,tsi) = E2D(mp.F4.compact.corr.inds);   % Actual field in estimation area
        Icube(:,:,tsi) = (abs(E2D).^2)*mp.WttlamVec(tsi)/mp.Wsum;
    end
    Isum2D = sum(Icube,3);


end %--END OF FUNCTION
