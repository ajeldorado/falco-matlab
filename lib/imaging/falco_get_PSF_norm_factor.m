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
%     figure(700+si); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im_temp_compact)); axis xy equal tight; colorbar; title('Compact Model'); drawnow; pause(0.2);
end

% %--Variables for DEBUGGING
% Ec = Etemp; 
% Ic = Im_temp_compact;


%--Compact Evaluation Model Normalizations
for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    Etemp = model_compact(mp, modvar,'NormOff','eval');
    Im_temp_compact = abs(Etemp).^2;
    mp.Fend.eval.I00(si) = max(max(Im_temp_compact));
%     figure(900+si); imagesc(log10(Im_temp_compact)); axis xy equal tight; colorbar; title('Compact EValuation Model'); drawnow; pause(0.2);
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
%         figure(800+counter); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Im_temp_full)); axis xy equal tight; colorbar; title('Full Model'); drawnow; pause(0.2);
        counter = counter+1;
    end
end

% %--Write files for DEBUGGING
% fitswrite(Ic/max(Ic(:)),'/Users/ajriggs/Downloads/compact_intensity.fits');
% fitswrite(real(Ec)/sqrt(max(Ic(:))),'/Users/ajriggs/Downloads/compact_real.fits');
% fitswrite(imag(Ec)/sqrt(max(Ic(:))),'/Users/ajriggs/Downloads/compact_imag.fits');
% fitswrite(Im_temp_full/max(Im_temp_full(:)),'/Users/ajriggs/Downloads/full_intensity.fits');
% fitswrite(real(Etemp)/sqrt(max(Im_temp_full(:))),'/Users/ajriggs/Downloads/full_real.fits');
% fitswrite(imag(Etemp)/sqrt(max(Im_temp_full(:))),'/Users/ajriggs/Downloads/full_imag.fits');
% keyboard

    
end %--END OF FUNCTION
