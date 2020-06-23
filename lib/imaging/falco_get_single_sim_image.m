% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get a simulated image from the full model at the specified
% wavelength and polarization state.
%
% ---------------
% INPUTS:
% - ic = index in list of combinations
% - inds_list = list of indices for the combinations
% - mp = structure of model parameters
%
% OUTPUTS
% - Iout: simulated image from full model at one wavelength and polarization
%
% REVISION HISTORY
% - Created on 2019-05-06 by A.J. Riggs.

function [Iout,varargout] = falco_get_single_sim_image(ic,inds_list,mp)

%--Get the starlight image
modvar.sbpIndex   = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(inds_list(1,ic)),1);
modvar.wpsbpIndex = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(inds_list(1,ic)),2);
mp.full.polaxis = mp.full.pol_conds(inds_list(2,ic));
modvar.whichSource = 'star';
if mp.fiber
    [Estar,Efiber] = model_full(mp, modvar);
    Ifiber = (abs(Efiber).^2);
    varargout{1} = Ifiber;
else
    Estar = model_full(mp, modvar);
end
Iout = (abs(Estar).^2); %--Apply spectral weighting outside this function
end %--END OF FUNCTION