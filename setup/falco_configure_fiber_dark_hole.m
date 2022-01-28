% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% FPM coordinates, [meters] and [dimensionless]

function mp = falco_configure_fiber_dark_hole(mp)

    V(1,1,:) = 2*pi*mp.fiber.a_phys*mp.fiber.NA./mp.sbp_centers;
    W = 1.1428*V - 0.996;
    U = sqrt(V.^2 - W.^2);
    
    maskFiberCore.rhoInner = 0;
    maskFiberCore.rhoOuter = mp.fiber.a;
    maskFiberCore.angDeg = 180;
    maskFiberCore.whichSide = mp.Fend.sides;
    
    maskFiberCladding.rhoInner = mp.fiber.a;
    maskFiberCladding.rhoOuter = 20;
    maskFiberCladding.angDeg = 180;
    maskFiberCladding.whichSide = mp.Fend.sides;
    
    if(mp.flagLenslet)
        %--Mask defining the area covered by the lenslet.  Only the immediate area
        %around the lenslet is propagated, saving computation time.  This lenslet
        %can then be moved around to different positions in Fend.
        maskLenslet.pixresFP = mp.Fend.res;
        maskLenslet.rhoInner = 0;
        maskLenslet.rhoOuter = mp.Fend.lensletWavRad;
        maskLenslet.angDeg = mp.Fend.corr.ang;
        maskLenslet.centering = mp.centering;
        maskLenslet.FOV = mp.Fend.FOV;
        maskLenslet.whichSide = mp.Fend.sides;
        [mp.Fend.lenslet.mask, ~, ~] = falco_gen_SW_mask(maskLenslet);
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
        %--Dummy mask to calculate the F5 coordinates correctly.
        maskF5.pixresFP = mp.F5.res;
        maskF5.rhoInner = 0;
        maskF5.rhoOuter = 1.22;
        maskF5.angDeg = 180;
        maskF5.centering = mp.centering;
        maskF5.FOV = mp.F5.FOV;
        maskF5.whichSide = mp.Fend.sides;
        [mp.F5.mask, mp.F5.xisDL, mp.F5.etasDL] = falco_gen_SW_mask(maskF5);
        
        %--Size of the output image in F5
        mp.F5.Nxi = size(mp.F5.mask, 2);
        mp.F5.Neta = size(mp.F5.mask, 1);
        
        %% Set up the fiber mode in F5
        
        maskFiberCore.pixresFP = mp.F5.res;
        maskFiberCore.FOV = mp.F5.FOV;
        [mp.F5.fiberCore.mask, ~, ~] = falco_gen_SW_mask(maskFiberCore);
        
        maskFiberCladding.pixresFP = mp.F5.res;
        maskFiberCladding.FOV = mp.F5.FOV;
        [mp.F5.fiberCladding.mask, ~, ~] = falco_gen_SW_mask(maskFiberCladding);
        
        [F5XIS, F5ETAS] = meshgrid(mp.F5.xisDL, mp.F5.etasDL);
        
        mp.F5.RHOS = sqrt((F5XIS - mp.F5.fiberPos(1)).^2 + (F5ETAS - mp.F5.fiberPos(2)).^2);
        mp.F5.fiberCore.mode = mp.F5.fiberCore.mask.*besselj(0, U.*mp.F5.RHOS/mp.fiber.a)./besselj(0,U);
        mp.F5.fiberCladding.mode = mp.F5.fiberCladding.mask.*besselk(0, W.*mp.F5.RHOS/mp.fiber.a)./besselk(0,W);
        mp.F5.fiberCladding.mode(isnan(mp.F5.fiberCladding.mode)) = 0;
        mp.F5.fiberMode = mp.F5.fiberCore.mode + mp.F5.fiberCladding.mode;
        fiberModeNorm = sqrt(sum(sum(abs(mp.F5.fiberMode).^2)));
        mp.F5.fiberMode = mp.F5.fiberMode./fiberModeNorm;
    else
                
        for nfib = 1:mp.Fend.Nfiber    
            maskFiberCore.pixresFP = mp.Fend.res;
            maskFiberCore.FOV = mp.Fend.FOV;
            maskFiberCore.xi_cen = mp.Fend.x_fiber(nfib);
            maskFiberCore.eta_cen = mp.Fend.y_fiber(nfib);
            [mp.Fend.fiberCore.mask, ~, ~] = falco_gen_SW_mask(maskFiberCore);

            maskFiberCladding.pixresFP = mp.Fend.res;
            maskFiberCladding.FOV = mp.Fend.FOV;
            maskFiberCladding.xi_cen = mp.Fend.x_fiber(nfib);
            maskFiberCladding.eta_cen = mp.Fend.y_fiber(nfib);
            [mp.Fend.fiberCladding.mask, ~, ~] = falco_gen_SW_mask(maskFiberCladding);

            [FENDXIS, FENDETAS] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
            
            FENDRHOS = sqrt((FENDXIS - mp.Fend.x_fiber(nfib)).^2 + (FENDETAS - mp.Fend.y_fiber(nfib)).^2);
            mp.Fend.fiberCore.mode = mp.Fend.fiberCore.mask.*besselj(0, U.*FENDRHOS/mp.fiber.a)./besselj(0,U);
            mp.Fend.fiberCladding.mode = mp.Fend.fiberCladding.mask.*besselk(0, W.*FENDRHOS/mp.fiber.a)./besselk(0,W);
            mp.Fend.fiberCladding.mode(isnan(mp.Fend.fiberCladding.mode)) = 0;
            fiberModeTemp = mp.Fend.fiberCore.mode + mp.Fend.fiberCladding.mode;
            fiberModeNorm = sqrt(sum(sum(abs(fiberModeTemp).^2)));
            mp.Fend.fiberMode(:,:,:,nfib) = fiberModeTemp./fiberModeNorm;
        end
        
    end

end %--END OF FUNCTION