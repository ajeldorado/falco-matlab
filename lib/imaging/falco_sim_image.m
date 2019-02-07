% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a simulated, single-sub-bandpass image with optional planet light added.
% -Measurement noise terms not added in yet
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


function Iout = falco_sim_image(mp, modvar, DM)


% facContrastToCounts = model_params.texp*model_params.peakCountsPerPixPerSec;
modvar.flagCalcJac = 0; % False for all cases
modvar.flagGetNormVal = false; % False for all cases
% % sbpIndex is already defined in modvar

if( strcmpi(modvar.whichSource,'offaxis') ) %--For throughput calculations by moving star off axis.
    modvar.wpsbpIndex = mp.wi_ref; %--At center wavelength only
    Eout = model_full(mp, DM, modvar);
    Iout = abs(Eout).^2;
    
else %--Regular
    Isum = 0; % Initialize unstacked image
    
%     for tsi=1:mp.jac.Nmode  %--Leave the tip/tilt part outside this function
        for wi=1:mp.Nwpsbp  % Add intensities from all incoherent sources separately
            % Starlight
            modvar.wpsbpIndex = wi;
            modvar.whichSource = 'star';
            Eout = model_full(mp, DM, modvar);
            Isum = Isum + (abs(Eout).^2);%*mp.jac.weights(tsi)/mp.Wsum;
        end 
%     end

    % Exoplanet light
    if(mp.planetFlag)
        for wi=1:mp.Nwpsbp % Calculate planet image and noise
                modvar.whichSource = 'exoplanet';
                modvar.lambdaIndex = wi;
                Eout = model_full(mp, DM,modvar);
                ImPlanetC = abs(Eout).^2; % In contrast
                Isum = Isum + ImPlanetC;
        end
    end

    Iout = Isum/mp.Nwpsbp;
end


%--OLD CODE: INCLUDES MEASUREMENT NOISE
% facContrastToCounts = model_params.texp*model_params.peakCountsPerPixPerSec;
% model_var.flagCalcJac = 0; % False for all cases
% model_var.flagGetNormVal = false; % False for all cases
% % % sbpIndex is already defined in model_var
% ImTot = 0; % Initialize unstacked image
% 
% % Stack num_im images:
% Isum = 0; % Initialize stacked image
% for qq = 1:model_params.num_im;
% 
% % % ImTot = 0; % Initialize unstacked image
% ImContrast = zeros(model_params.NetaCam,model_params.NxiCam,model_params.Nwpsbp);
% for wi=1:model_params.Nwpsbp  % Add intensities from all incoherent sources separately
%     % Starlight
%     model_var.wpsbpIndex = wi;
%     model_var.whichSource = 'star';
%     Eout = model_full(model_params, DM_config, model_var);
%     ImContrast(:,:,wi) = abs(Eout).^2; % In units of contrast
% end    
% 
% % % Convert to counts:
% % DHmeanCountsBefore = sum(sum( mean(ImContrast,3).*model_params.ScoreMask))/sum(sum(model_params.ScoreMask));
% 
% 
% for wi=1:model_params.Nwpsbp % Add photon noise to starlight
% 
%     ImMasked = ImContrast(:,:,wi).*model_params.fieldStop;
%     ImCounts = ImMasked*facContrastToCounts/model_params.Nwpsbp; 
%     % Add photon noise here
%     if(model_params.noiseFlag)
%         ImCountsShotNoise = round( 1e12*imnoise(ImCounts*1e-12,'poisson') );% round after adding photon noise
%     else
%         ImCountsShotNoise = ImCounts;
%     end
%     ImTot = ImTot+ImCountsShotNoise;
% end
% 
% % Planet light
% if(model_params.planetFlag)
%     for wi=1:model_params.Nwpsbp % Calculate planet image and noise
% 
%             model_var.whichSource = 'exoplanet';
%             model_var.lambdaIndex = wi;
%             Eout = model_full(model_params, DM_config,model_var);
%             ImC = abs(Eout).^2; % In contrast
% 
%             % Convert to counts:
%             ImMasked = ImC.*model_params.fieldStop;
%             ImCounts = (ImMasked*facContrastToCounts)/model_params.Nwpsbp; 
%             % Add photon noise here
%             if(model_params.noiseFlag)
%                 ImCountsShotNoise = round( 1e12*imnoise(ImCounts*1e-12,'poisson') );% round after adding photon noise
%                 % Should the rounding occur before or after the noisy planet
%                 % light is added to the starlight?
%             else
%                 ImCountsShotNoise = ImCounts;
%             end
%             ImTot = ImTot+ImCountsShotNoise;
%     end
% end
%     
% % Add zodi noise here
% if(model_params.zodiFlag)
%     ImZodiCounts = model_params.c_zodi*model_params.fieldStop*facContrastToCounts;
%     if(model_params.noiseFlag)
%         ImZodiCountsShotNoise = round( 1e12*imnoise(ImZodiCounts*1e-12,'poisson') );
%     else
%         ImZodiCountsShotNoise = ImZodiCounts;
%     end
%     ImTot = ImTot + ImZodiCountsShotNoise;
% end
% % ImTot = ImTot/model_params.num_im/model_params.Nwpsbp; %--Need to do this
% % before adding the Poisson noise!
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Add readout noise here
% if(model_params.noiseFlag && model_params.readNoiseStd>0)
%     meanVal = 0;
%     %readNoise = (round(random('Normal',meanVal,model_params.readNoiseStd,[model_params.NetaCam, model_params.NxiCam]))); %--random.m requires the stats toolbox
%     readNoise = round(model_params.readNoiseStd*randn(model_params.NetaCam, model_params.NxiCam));
% else
%     readNoise = 0;
% end
% 
% 
% ImNoisy = ImTot+readNoise;  % Add readout noise in terms of counts
% Isum = Isum+ImNoisy/facContrastToCounts;  % Convert back to contrast
% end  % End of image stacking
% 
% Iout = Isum/model_params.num_im;
% 
