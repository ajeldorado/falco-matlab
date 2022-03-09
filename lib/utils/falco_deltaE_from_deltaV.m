% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute the model-based delta E-field for the given delta DM voltages at
% the specified wavelength.
%
% INPUTS
% ------
% mp : structure of model parameters
% dV1 : 2-D map of delta commands to apply to DM1
% dV2 : 2-D map of delta commands to apply to DM2
% lambda : wavelength at which to evaluate the PSF
%
% OUTPUTS
% -------
% deltaE : normalized 2-D complex delta E-field at final focus
 
function deltaE = falco_deltaE_from_deltaV(mp, dV1, dV2, lambda)
 
    % Store initial DM settings
    V1init = mp.dm1.V;
    V2init = mp.dm2.V;
 
    modvar.lambda = lambda;
    modvar.whichSource = 'star';
    modvar.starIndex = 1; % ignored because of modvar.lambda
    modvar.zernIndex = 1; % ignored because of modvar.lambda
    modvar.sbpIndex = 1;
    modvar.wpsbpIndex = 1;
 
    EforNorm = model_compact(mp, modvar, 'getnorm');
    I00 = max(max(abs(EforNorm).^2));
 
    Epre = model_compact(mp, modvar, 'nonorm') / sqrt(I00);
 
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, V1init+dV1);
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, V2init+dV2);
    Epost = model_compact(mp, modvar, 'nonorm') / sqrt(I00);
    deltaE = Epost - Epre;
 
    % Reset to original voltages
    mp.dm1.V = V1init;
    mp.dm2.V = V2init;
 
    if mp.flagPlot
        figure(151); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, abs(deltaE));
        axis xy equal tight; colorbar; colormap parula;
        title('abs(deltaE)', 'Fontsize', 20);
        xlabel('\lambda/D');
        ylabel('\lambda/D');
        set(gca, 'Fontsize', 20);
        drawnow;
 
        figure(152); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, angle(deltaE));
        axis xy equal tight; colorbar; colormap hsv;
        title('angle(deltaE)', 'Fontsize', 20)
        xlabel('\lambda/D');
        ylabel('\lambda/D');
        set(gca, 'Fontsize', 20);
        drawnow;
    end
 
end

