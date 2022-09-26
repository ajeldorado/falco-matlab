% Copyright 2018, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Measure the cross term of the E-field and the probe in a single Jacobian mode. Used for IEFC.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% iefcStruct : structure of output variables for IEFC

function iefcStruct = falco_est_delta_intensity(mp, whichDM, dV, iJacMode)

    flagPlot = false;  % for debugging
    
    modvar = ModelVariables;
    modvar.whichSource = 'star';
    modvar.sbpIndex = mp.jac.sbp_inds(iJacMode);
    modvar.starIndex = mp.jac.star_inds(iJacMode);
    modvar.zernIndex = mp.jac.zern_inds(iJacMode);
  
    if whichDM == 1
        
        V0 = mp.dm1.V;
        mp.dm1.V = V0 + dV;
        imagePlus = falco_get_sbp_image(mp, modvar.sbpIndex);
        mp.dm1.V = V0 - dV;
        imageMinus = falco_get_sbp_image(mp, modvar.sbpIndex);
        mp.dm1.V = V0; % reset;
        
    elseif whichDM == 2
        
        V0 = mp.dm2.V;
        mp.dm2.V = V0 + dV;
        imagePlus = falco_get_sbp_image(mp, modvar.sbpIndex);
        mp.dm2.V = V0 - dV;
        imageMinus = falco_get_sbp_image(mp, modvar.sbpIndex);
        mp.dm2.V = V0; % reset;
        
    end
    DI2D = (imagePlus - imageMinus);

    if flagPlot
        figure(1001); imagesc(log10(imagePlus)); axis xy equal tight; colorbar; drawnow;
        figure(1002); imagesc(log10(imageMinus)); axis xy equal tight; colorbar; drawnow;

        figure(1005); imagesc(log10(abs(DI2D))); axis xy equal tight; colorbar; drawnow;
        figure(1006); imagesc(angle(DI2D)); axis xy equal tight; colorbar; drawnow;
    end

    iefcStruct.DeltaI = DI2D(mp.Fend.corr.maskBool).';
    iefcStruct.imagePlus = imagePlus;
    iefcStruct.imageMinus = imageMinus;
    
end %--END OF FUNCTION
