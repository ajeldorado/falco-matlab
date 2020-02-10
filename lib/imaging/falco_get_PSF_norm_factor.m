% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to compute the normalization value for the compact and full
% models. 
% For vortex coronagraphs: The normalization value is the peak intensity 
%   for an on-axis object with the entire coronagraph in place except with 
%   the focal plane mask removed.
% For other coronagraphs: The normalization value is the peak intensity 
%   for an off-axis object at separation (mp.source_x_offset_norm,
%   mp.source_y_offset_norm) lambda0/D.
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
%
% OUTPUTS
% - mp = structure of model parameters

function mp = falco_get_PSF_norm_factor(mp)

%--Different normalization factor used when comparing to PROPER model:
mp.sumPupil = sum(sum(abs(mp.P1.compact.mask.*padOrCropEven(mean(mp.P1.compact.E,3),size(mp.P1.compact.mask,1) )).^2));

%--Initialize Model Normalizations
mp.Fend.compact.I00 = ones(1,mp.Nsbp); % Initial input before computing
mp.Fend.eval.I00 = ones(1,mp.Nsbp); % Initial input before computing
mp.Fend.full.I00 = ones(mp.Nsbp,mp.Nwpsbp); % Initial input before computing

modvar.zernIndex = 1;
modvar.whichSource = 'star';    

%--Compact Model Normalizations
for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    Etemp = model_compact(mp, modvar,'getNorm');
    mp.Fend.compact.I00(si) = max(max(abs(Etemp).^2));
end

%--Compact Evaluation Model Normalizations
for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    Etemp = model_compact(mp, modvar,'getNorm','eval');
    mp.Fend.eval.I00(si) = max(max(abs(Etemp).^2));
end

%--Full Model Normalizations (at points for entire-bandpass evaluation)
if(mp.flagSim)
    if(mp.flagParfor)
        parfor li = 1:mp.Nsbp*mp.Nwpsbp
            I00vec{li} = model_full_norm_wrapper(li,mp);
        end
        
        counter = 0;
        for si=1:mp.Nsbp
            for wi=1:mp.Nwpsbp
                counter = counter+1;
                mp.Fend.full.I00(si,wi) = I00vec{counter};
            end
        end
    else %--No parfor
        for si=1:mp.Nsbp
            for wi=1:mp.Nwpsbp
                modvar.sbpIndex = si;
                modvar.wpsbpIndex = wi;
                Etemp = model_full(mp, modvar,'getNorm');
                mp.Fend.full.I00(si,wi) = max(max(abs(Etemp).^2));
            end
        end
    end
    
end

%--Visually verify the normalized coronagraphic PSF
modvar.ttIndex = 1;
modvar.sbpIndex = mp.si_ref;
modvar.wpsbpIndex = mp.wi_ref;
modvar.whichSource = 'star';     
E0c = model_compact(mp, modvar);
I0c = abs(E0c).^2;
if(mp.flagPlot)
    figure(501); imagesc(log10(I0c)); axis xy equal tight; colorbar;
    title('(Compact Model: Normalization Check Using Starting PSF)'); 
    drawnow;
end
E0f = model_full(mp, modvar);
I0f = abs(E0f).^2;
if(mp.flagPlot)
    figure(502); imagesc(log10(I0f)); axis xy equal tight; colorbar;
    title('(Full Model: Normalization Check Using Starting PSF)'); drawnow;
end

end %--END OF FUNCTION


%--Extra function needed to use parfor (because parfor can have only a
%  single changing input argument).
function I00 = model_full_norm_wrapper(li,mp)
    modvar.sbpIndex = mp.full.indsLambdaMat(li,1);
    modvar.wpsbpIndex = mp.full.indsLambdaMat(li,2);
    modvar.zernIndex = 1;
    modvar.whichSource = 'star'; 
    
    Etemp = model_full(mp, modvar,'getNorm');
    I00 = max(max(abs(Etemp).^2));
end