% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Get a simulated image from the full model at the specified
% wavelength and polarization state and for the specified star.
%
% INPUTS
% ------
% ic : index in list of combinations
% inds_list : list of indices for the combinations
% mp : structure of model parameters
%
% OUTPUTS
% -------
% Iout : simulated image from full model at one wavelength and polarization
% and for a single specified star

function Iout = falco_get_single_sim_image(ic, inds_list, mp)

    modvar.sbpIndex   = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(inds_list(1, ic)), 1);
    modvar.wpsbpIndex = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(inds_list(1, ic)), 2);
    modvar.starIndex = inds_list(3, ic);
    mp.full.polaxis = mp.full.pol_conds(inds_list(2, ic));
    modvar.whichSource = 'star';
    Estar = model_full(mp, modvar);
    Iout = abs(Estar).^2; %--Apply spectral and stellar weighting outside this function
    
end
