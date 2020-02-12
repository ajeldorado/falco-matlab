% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Resolution at Final Focal Plane (Fend)

function mp = falco_get_Fend_resolution(mp)

    mp.Fend.dxi = (mp.fl*mp.lambda0/mp.P4.D)/mp.Fend.res; % sampling at Fend.[meters]
    mp.Fend.deta = mp.Fend.dxi; % sampling at Fend.[meters]    
    if(mp.flagLenslet)
        mp.Fend.lenslet.D = 2*mp.Fend.res*mp.Fend.lensletWavRad*mp.Fend.dxi;
        mp.Fend.x_lenslet_phys = mp.Fend.dxi*mp.Fend.res*mp.Fend.x_lenslet;
        mp.Fend.y_lenslet_phys = mp.Fend.deta*mp.Fend.res*mp.Fend.y_lenslet;

        mp.F5.dxi = mp.lensletFL*mp.lambda0/mp.Fend.lenslet.D/mp.F5.res;
        mp.F5.deta = mp.F5.dxi;
    end

    %--Compact evaluation model at higher resolution
    mp.Fend.eval.dxi = (mp.fl*mp.lambda0/mp.P4.D)/mp.Fend.eval.res; % [meters]
    mp.Fend.eval.deta = mp.Fend.eval.dxi; % [meters]  
    
end %--END OF FUNCTION