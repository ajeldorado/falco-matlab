% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to compute the normalization value for the compact and full
% models. The normalization value is the peak intensity for an on-axis
% object with the entire coronagraph in place except with the focal plane
% mask removed.
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
%
% OUTPUTS
% - mp = structure of model parameters
%
% REVISION HISTORY: 
% -Created on 2018-01-24 by A.J. Riggs.


function mp = falco_get_PSF_norm_factor(mp)

%--Initialize Model Normalizations
mp.Fend.compact.I00 = ones(1,mp.Nsbp); % Initial input before computing
mp.Fend.eval.I00 = ones(1,mp.Nsbp); % Initial input before computing
mp.Fend.full.I00 = ones(1,mp.full.Nlam); % Initial input before computing

modvar.zernIndex = 1;
% modvar.wpsbpIndex = wi;
% modvar.ttIndex = 1; % 1 is the zero-offset tip/tilt setting
modvar.whichSource = 'star';    

%--Compact Model Normalizations
for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    Etemp = model_compact(mp, modvar,'NormOff');
    Im_temp_compact = abs(Etemp).^2;
    mp.Fend.compact.I00(si) = max(max(Im_temp_compact));
    %     figure; imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(mean(Im_temp_compact,3)/mp.Fend.compact.I00(si)),[-9 0]); axis xy equal tight; colorbar; title('Compact Model');
end

%--Compact Evaluation Model Normalizations
for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    Etemp = model_compact(mp, modvar,'NormOff','eval');
    Im_temp_compact = abs(Etemp).^2;
    mp.Fend.eval.I00(si) = max(max(Im_temp_compact));
    %     figure; imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(mean(Im_temp_compact,3)/mp.Fend.compact.I00(si)),[-9 0]); axis xy equal tight; colorbar; title('Compact Model');
end

% %--Full Model Normalizations (at points for entire-bandpass evaluation)
% for ilam=1:mp.full.Nlam
%     %modvar.lambda = mp.full.lambdas(ilam);
%     Etemp = model_full(mp, modvar,'NormOff');
%     %Etemp = model_compact(mp, modvar,'NormOff');
%     Im_temp_full = abs(Etemp).^2;
%     mp.Fend.full.I00(ilam) = max(max(Im_temp_full));
% end

%--Full Model Normalizations (at points for entire-bandpass evaluation)
counter = 1;
for si=1:mp.Nsbp
    for wi=1:mp.Nwpsbp
        modvar.sbpIndex = si;
        modvar.wpsbpIndex = wi;
        Etemp = model_full(mp, modvar,'NormOff');
        Im_temp_full = abs(Etemp).^2;
        mp.Fend.full.I00(counter) = max(Im_temp_full(:));
        counter = counter+1;
    end
end

    
end %--END OF FUNCTION
