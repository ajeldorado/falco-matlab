% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Crop off extra zero padding around the Lyot stop to speed up MFT propagation.
% For both the full and compact models.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_crop_lyot_stop(mp)

    % Full model
    if ~mp.full.flagPROPER
        LSsum = sum(mp.P4.full.mask(:));
        LSdiff = 0;
        counter = 2;
        while(abs(LSdiff) <= 1e-7)
            mp.P4.full.Narr = length(mp.P4.full.mask) - counter; % number of points across the cropped-down Lyot stop
            LSdiff = LSsum - sum(sum(padOrCropEven(mp.P4.full.mask, mp.P4.full.Narr-2))); %--Subtract an extra 2 to negate the extra step that overshoots.
            counter = counter + 2;
        end
        mp.P4.full.croppedMask = padOrCropEven(mp.P4.full.mask, mp.P4.full.Narr); %--The cropped-down Lyot stop for the full model. 
    end

    % Compact model
    LSsum = sum(mp.P4.compact.mask(:));
    LSdiff = 0;
    counter = 2;
    while(abs(LSdiff) <= 1e-7)
        mp.P4.compact.Narr = length(mp.P4.compact.mask) - counter; % number of points across the cropped-down Lyot stop
        LSdiff = LSsum - sum(sum(padOrCropEven(mp.P4.compact.mask, mp.P4.compact.Narr-2))); %--Subtract an extra 2 to negate the extra step that overshoots.
        counter = counter + 2;
    end
    mp.P4.compact.croppedMask = padOrCropEven(mp.P4.compact.mask, mp.P4.compact.Narr); %--The cropped-down Lyot stop for the compact model

end
