% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass.
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - si = index of sub-bandpass for which to take the image
%
% OUTPUTS
% - Isum: sub-bandpass image
%
% REVISION HISTORY
% - Created on 2019-02-06 by A.J. Riggs.


function Isum = falco_get_sim_sbp_image(mp,si)


%--Compute the DM surfaces outside the full model to save lots of time
if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end

Isum = 0; % Initialize image
modvar.sbpIndex = si;
% facContrastToCounts = model_params.texp*model_params.peakCountsPerPixPerSec;

%--Get the starlight image
modvar.whichSource = 'star';
for wi=1:mp.Nwpsbp
    modvar.wpsbpIndex = wi;
    Etemp = model_full(mp, modvar);
    Isum = Isum + (abs(Etemp).^2)*mp.full.lambda_weights(wi); %--Do not apply sbp_weight unless full bandpass image is being created.
end 

%--Include the planet image if flagged
if(mp.planetFlag)
    modvar.whichSource = 'exoplanet';
    for wi=1:mp.Nwpsbp
        modvar.wpsbpIndex = wi;
        Eplanet = model_full(mp,modvar);
        Isum = Isum + abs(Eplanet).^2*mp.full.lambda_weights(wi); %--Do not apply sbp_weight unless full bandpass image is being created.
    end
end
    
    
end %--END OF FUNCTION

