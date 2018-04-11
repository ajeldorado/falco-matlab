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
% REVISION HISTORY: 
% -Created on 2018-01-24 by A.J. Riggs.


function mp = falco_get_PSF_norm_factor(mp, DM)

%--Initialize Model Normalizations
mp.F4.compact.I00 = ones(1,mp.Nsbp); % Initial input before computing
mp.F4.full.I00 = ones(1,mp.Nsbp); % Initial input before computing

%--Compute Model Normalizations
modvar.flagGetNormVal=true; %--Turn off normalization mode in compact and full models
for si=1:mp.Nsbp,
    Im_temp_full = zeros(mp.F4.full.Neta,mp.F4.full.Nxi,mp.Nwpsbp);

    for wi=1:mp.Nwpsbp
        modvar.flagCalcJac = 0; 
        modvar.sbpIndex = si;
        modvar.wpsbpIndex = wi;
        modvar.ttIndex = 1; % 1 is the zero-offset tip/tilt setting
        modvar.whichSource = 'star';     
        
        EtempFull = model_full(mp, DM, modvar);
        Im_temp_full(:,:,wi) = abs(EtempFull).^2;
    end
    
    EtempCompact = model_compact(mp, DM, modvar);
    Im_temp_compact = abs(EtempCompact).^2;
        
    mp.F4.full.I00(si) = max(max(mean(Im_temp_full,3)));
    mp.F4.compact.I00(si) = max(max(Im_temp_compact));
%     figure; imagesc(mp.F4.full.xisDL,mp.F4.full.etasDL,log10(mean(Im_temp_full,3)/mp.F4.full.I00(si)),[-9 0]); axis xy equal tight; colorbar; title('Full Model');
%     figure; imagesc(mp.F4.compact.xisDL,mp.F4.compact.etasDL,log10(mean(Im_temp_compact,3)/mp.F4.compact.I00(si)),[-9 0]); axis xy equal tight; colorbar; title('Compact Model');
end
modvar.flagGetNormVal=false; %--Turn off normalization mode in the models

end
