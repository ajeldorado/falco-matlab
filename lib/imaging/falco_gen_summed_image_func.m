% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a simulated broadband image with the full model at
% the weighted wavelengths specified in the structure mp.entireBWimage. No
% tip/tilt offsets are included.
%
% X Function to generate a simulated broadband image the specified 
%  wavelength-tip/tilt setting given (by "tsi").
%  -> To be called by a wrapper function. Set up as a function in order to
%  allow use of parfor.
%
% Modified on 2018-09-17 by A.J. Riggs to remove mp.entireBWimage and to
% use the normalized weights.
% Modified on 2018-08-31 by A.J. Riggs to use the structure mp.entireBWimage and not to include tip/tilt
% Modified on 2018-05-17 by A.J. Riggs to be a function and be for just one
% wavelength-tip/tilt pair instead of all.
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


function Isum = falco_gen_summed_image_func(si,mp, DM)


    % facContrastToCounts = model_params.texp*model_params.peakCountsPerPixPerSec;
    modvar.flagCalcJac = 0; % False for all cases
    modvar.flagGetNormVal = false; % False for all cases
    modvar.sbpIndex = si;%mp.Wttlam_si(tsi);
%     modvar.ttIndex = mp.Wttlam_ti(tsi);

    Isum = 0; % Initialize image

    %--Loop over wavelengths within each sub-bandpass
    for wi=1:mp.Nwpsbp  % Add intensities from all incoherent sources separately
        modvar.wpsbpIndex = wi;

        %--Starlight
        modvar.whichSource = 'star';
        Eout = model_full(mp, DM, modvar);
        Isum = Isum + (abs(Eout).^2)*mp.full.lambda_weights(wi);

        %--Exoplanet light
        if(mp.planetFlag)
            modvar.whichSource = 'exoplanet';
            Eout = model_full(mp, DM,modvar);
            ImPlanetC = abs(Eout).^2; % In contrast
            Isum = Isum + ImPlanetC*mp.full.lambda_weights(wi);
        end

    end 

    
% %     %--Another set of wavelengths only for use in generating the image of the
% %     % entire-bandwdith image with the full model. 
% %     mp.full.Nlam = mp.Nsbp;
% %     mp.full.lambdas = mp.lambda0*linspace( 1-mp.fracBW/2,1+mp.fracBW/2,mp.full.Nlam);
% %     mp.entireBWimage.Iweights = ones(mp.full.Nlam,1);
% %     if(mp.full.Nlam>2)
% %         mp.entireBWimage.Iweights(1) = 1/2;
% %         mp.entireBWimage.Iweights(end) = 1/2;
% %     end
    
% %     IweightsSum = sum(mp.entireBWimage.Iweights);
% % 
% %     %--Loop over wavelengths within each sub-bandpass
% %     for ilam = 1:mp.full.Nlam
% % 
% %         %--Starlight
% %         modvar.lambda = mp.full.lambdas(ilam);
% %         modvar.whichSource = 'star';
% %         Eout = model_full(mp, DM, modvar);
% %         Isum = Isum + mp.entireBWimage.Iweights(ilam)/IweightsSum*(abs(Eout).^2);
% % %         Isum = Isum + (abs(Eout).^2)*(mp.jac.weights(tsi)/mp.Wsum)/mp.Nwpsbp;
% %         
% % %         %--Exoplanet light
% % %         if(mp.planetFlag)
% % %             modvar.whichSource = 'exoplanet';
% % %             Eout = model_full(mp, DM,modvar);
% % %             ImPlanetC = abs(Eout).^2; % In contrast
% % %             Isum = Isum + ImPlanetC*(mp.jac.weights(tsi)/mp.Wsum)/mp.Nwpsbp;
% % %         end
% %     end
    
    

    
end %--END OF FUNCTION

