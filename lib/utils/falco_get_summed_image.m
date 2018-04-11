% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a simulated broadband image for all wavelengths and
% tip/tilt settings provided.
%  -> Purpose is to be easier to use and faster than re-computing DM
%  surfaces for each call of model_full.m
%
% Modified on 23 May 2017 by A.J. Riggs to not include any noise sources.
%   Now intended only for coronagraph design. New function name is
%   falco_get_image.m.
% Changed on 19 September 2016 by A.J. Riggs from hcil_* to
%  fpwc_getSimImage.m
% 7 May 2015: Fixed error in which I was zeroing out ImTot instead of
%  stacking it. Moved ImTot=0 outside qq loop and inserted line to divide
%  ImTot by number of images at the end.
%  getSimImage_HCIL_BB_v1.m - 26 Nov. 2014: 
%   -Now uses a flux and exposure time instead of desired avg counts to set
%   number of counts.
% getSimImage_HCIT_BB_v2.m - 23 Sept. 2014: Modified to create an image from several discrete
%  wavelengths instead of just 1. Also added a separate flag for the
%  planet.


function Isum = falco_get_summed_image(mp, DM)

%--Compute the DM surfaces outside the full model to save lots of time
if(any(DM.dm_ind==1)); DM.dm1.surfM = falco_gen_dm_surf(DM.dm1,DM.dm1.dx,DM.dm1.NdmPad); end;
if(any(DM.dm_ind==2)); DM.dm2.surfM = falco_gen_dm_surf(DM.dm2,DM.dm2.dx,DM.dm2.NdmPad); end;

% facContrastToCounts = model_params.texp*model_params.peakCountsPerPixPerSec;
modvar.flagCalcJac = 0; % False for all cases
modvar.flagGetNormVal = false; % False for all cases
     
Isum = 0; % Initialize image

%--Loop over pairs of sub-bandpass and tip/tilt offsets
for tsi=1:mp.Nttlam
    modvar.sbpIndex = mp.Wttlam_si(tsi);
    modvar.ttIndex = mp.Wttlam_ti(tsi);
    
    %--Loop over wavelengths within each sub-bandpass
    for wi=1:mp.Nwpsbp  % Add intensities from all incoherent sources separately
        modvar.wpsbpIndex = wi;

        %--Starlight
        modvar.whichSource = 'star';
        Eout = model_full(mp, DM, modvar);
        Isum = Isum + (abs(Eout).^2)*(mp.WttlamVec(tsi)/mp.Wsum)/mp.Nwpsbp;

        %--Exoplanet light
        if(mp.planetFlag)
            modvar.whichSource = 'exoplanet';
            Eout = model_full(mp, DM,modvar);
            ImPlanetC = abs(Eout).^2; % In contrast
            Isum = Isum + ImPlanetC*(mp.WttlamVec(tsi)/mp.Wsum)/mp.Nwpsbp;
        end

    end 
end
    
end %--END OF FUNCTION

