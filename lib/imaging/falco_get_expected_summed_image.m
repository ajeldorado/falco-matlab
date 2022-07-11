% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Ibandavg = falco_get_expected_summed_image(mp,cvar)
%
% Function to generate the expected broadband image over the entire
% bandpass by adding the model-based delta electric field on top of the
% current E-field estimate in each sub-bandpass.
%
% INPUTS
% mp : structure of all model parameters
%
% OUTPUTS
% Ibandavg : band-averaged image in units of normalized intensity

function Ibandavg = falco_get_expected_summed_image(mp, cvar)

    %--Initialize variables
    Ibandavg = 0;
    EnewTempVecArray = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
    EoldTempVecArray = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);

    %--Generate the model-based E-field with the new DM setting
    modvar = ModelVariables;
    for iMode = 1:mp.jac.Nmode    
        modvar.sbpIndex = mp.jac.sbp_inds(iMode);
        modvar.zernIndex = mp.jac.zern_inds(iMode);
        modvar.starIndex = mp.jac.star_inds(iMode);
        modvar.whichSource = 'star';
        modvar.wpsbpIndex = 0; %--Dummy, placeholder value
        Etemp = model_compact(mp, modvar);
        EnewTempVecArray(:, iMode) = Etemp(mp.Fend.corr.mask);
    end
    
    %--Revert to the previous DM commands
    if(any(mp.dm_ind==1));  mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V - mp.dm1.dV);  end
    if(any(mp.dm_ind==2));  mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V - mp.dm2.dV);  end
    if(any(mp.dm_ind==3));  mp.dm3.V = mp.dm3.V - mp.dm3.dV;  end
    if(any(mp.dm_ind==4));  mp.dm4.V = mp.dm4.V - mp.dm4.dV;  end
    if(any(mp.dm_ind==5));  mp.dm5.V = mp.dm5.V - mp.dm5.dV;  end
    if(any(mp.dm_ind==6));  mp.dm6.V = mp.dm6.V - mp.dm6.dV;  end
    if(any(mp.dm_ind==7));  mp.dm7.V = mp.dm7.V - mp.dm7.dV;  end
    if(any(mp.dm_ind==8));  mp.dm8.V = mp.dm8.V - mp.dm8.dV;  end
    if(any(mp.dm_ind==9));  mp.dm9.V = mp.dm9.V - mp.dm9.dV;  end    
        
    %--Generate the model-based E-field with the previous DM setting
    for iMode = 1:mp.jac.Nmode    
        modvar.sbpIndex = mp.jac.sbp_inds(iMode);
        modvar.zernIndex = mp.jac.zern_inds(iMode);
        modvar.starIndex = mp.jac.star_inds(iMode);
        modvar.whichSource = 'star';
        modvar.wpsbpIndex = 0; %--Dummy, placeholder value
        Etemp = model_compact(mp, modvar);
        EoldTempVecArray(:, iMode) = Etemp(mp.Fend.corr.mask);
    end
    
    %--Compute the expected new 2-D intensity image
    for iMode = 1:mp.jac.Nmode    
        EexpectedVec = cvar.Eest(:, iMode) + (EnewTempVecArray(:, iMode) - EoldTempVecArray(:, iMode));
        Eexpected2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
        Eexpected2D(mp.Fend.corr.mask) = EexpectedVec;
        
        sbpIndex = mp.jac.sbp_inds(iMode);
        Ibandavg = Ibandavg +  mp.sbp_weights(sbpIndex)*abs(Eexpected2D).^2;
    end
    
end %--END OF FUNCTION
