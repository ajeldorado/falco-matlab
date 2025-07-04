% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute the normalization value for the compact and full models.
%
% For vortex coronagraphs: The normalization value is the peak intensity 
%   for an on-axis object with the entire coronagraph in place except with 
%   the focal plane mask removed.
% For other coronagraphs: The normalization value is the peak intensity 
%   for an off-axis object at separation (mp.source_x_offset_norm,
%   mp.source_y_offset_norm) lambda0/D.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_compute_psf_norm_factor(mp, varargin)

%--If there is an extra input, it is for zeroing out the (delta) DM voltages temporarily.
zeroOutVoltages = false; % default
if size(varargin, 2) == 1
    extraArg = varargin{1};
    if extraArg == 0
        zeroOutVoltages = true;
    else
        error('The 2nd argument must be 0 to zero out the DM voltages.')
    end
end

if zeroOutVoltages
    V1init = mp.dm1.V;
    V2init = mp.dm2.V;
    mp.dm1.V = 0*mp.dm1.V;
    mp.dm2.V = 0*mp.dm2.V;
end


%--Different normalization factor used when comparing to PROPER model:
mp.sumPupil = sum(sum(abs(mp.P1.compact.mask.*padOrCropEven(mean(mp.P1.compact.E,3),size(mp.P1.compact.mask,1) )).^2));

%--Initialize Model Normalizations
mp.Fend.compact.I00 = ones(1, mp.Nsbp); % Initial input before computing
mp.Fend.eval.I00 = ones(1, mp.Nsbp); % Initial input before computing
mp.Fend.full.I00 = ones(mp.Nsbp, mp.Nwpsbp); % Initial input before computing
if mp.flagFiber
    mp.Fend.compact.I00_fiber = ones(mp.Fend.Nfiber, mp.Nsbp); % Initial input before computing
    mp.Fend.eval.I00_fiber = ones(mp.Fend.Nfiber, mp.Nsbp); % Initial input before computing
    mp.Fend.full.I00_fiber = ones(mp.Fend.Nfiber,mp.Nsbp, mp.Nwpsbp); % Initial input before computing
end

modvar = ModelVariables;
modvar.zernIndex = 1;
modvar.whichSource = 'star';  
modvar.starIndex = 1; % Always use first star for image normalization

%--Compact Model Normalizations
for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    if mp.flagFiber
        [Etemp,Efibertemp] = model_compact(mp, modvar, 'getNorm');
        mp.Fend.compact.I00_fiber(:,si) = abs(Efibertemp).^2;
    else
        Etemp = model_compact(mp, modvar, 'getNorm');
    end
    mp.Fend.compact.I00(si) = max(max(abs(Etemp).^2));
end

%--Compact Evaluation Model Normalizations
for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    if mp.flagFiber
%         [Etemp,Efibertemp] = model_compact(mp, modvar, 'getNorm', 'eval');
        mp.Fend.eval.I00_fiber(:,si) = mp.Fend.compact.I00_fiber(:,si);
    else
        Etemp = model_compact(mp, modvar, 'getNorm', 'eval');
    end
    mp.Fend.eval.I00(si) = max(max(abs(Etemp).^2));
end

%--Full Model Normalizations (at points for entire-bandpass evaluation)
if(mp.flagParfor)

    % Remove testbed object
    if isfield(mp, 'tb')
        mptmp = rmfield(mp, 'tb');
    else
        mptmp = mp;
    end
    
    parfor li = 1:mp.Nsbp*mp.Nwpsbp
        if mptmp.flagFiber
            [I00vec{li},I00fibervec{li}] = model_full_norm_wrapper(li,mptmp);
        else
            I00vec{li} = model_full_norm_wrapper(li,mptmp);
        end
    end

    counter = 0;
    for si=1:mp.Nsbp
        for wi=1:mp.Nwpsbp
            counter = counter+1;
            if mptmp.flagFiber
                mp.Fend.full.I00_fiber(si, wi) = I00fibervec{counter};
            end
            mp.Fend.full.I00(si, wi) = I00vec{counter};
        end
    end
    
else %--No parfor
    for si=1:mp.Nsbp
        for wi=1:mp.Nwpsbp
            modvar.sbpIndex = si;
            modvar.wpsbpIndex = wi;
            if mp.flagFiber
                [Etemp,Efibertemp] = model_full(mp, modvar, 'getNorm');
                mp.Fend.full.I00_fiber(:,si) = abs(Efibertemp).^2;
            else
                Etemp = model_full(mp, modvar,'getNorm');
            end
            mp.Fend.full.I00(si,wi) = max(max(abs(Etemp).^2));
        end
    end
end

%--Visually verify the normalized coronagraphic PSF
if(mp.flagPlot)
    modvar.sbpIndex = mp.si_ref;
    modvar.wpsbpIndex = mp.wi_ref;
    modvar.zernIndex = 1;
    modvar.starIndex = 1; % Always use first star for image normalization
    modvar.whichSource = 'star'; 

    E0c = model_compact(mp, modvar);
    I0c = abs(E0c).^2;
    figure(501); imagesc(log10(I0c)); axis xy equal tight; colorbar;
    title('Compact Model: Normalization Check');
    set(gca,'Fontsize', 18)
    drawnow;

    E0f = model_full(mp, modvar);
    I0f = abs(E0f).^2;
    figure(502); imagesc(log10(I0f)); axis xy equal tight; colorbar;
    title('Full Model: Normalization Check');
    set(gca,'Fontsize', 18)
    drawnow;
end

if zeroOutVoltages

    % Reset voltages to original value
    mp.dm1.V = V1init;
    mp.dm2.V = V2init;

    % Save these normalizations to other arrays
    mp.I00compactZeros = mp.Fend.compact.I00(:);
    mp.I00fullZeros = mp.Fend.full.I00;

else

    mp.Fend.compact.I00ratio = mp.Fend.compact.I00(:) ./ mp.I00compactZeros(:);
    mp.Fend.full.I00ratio = mp.Fend.full.I00 ./ mp.I00fullZeros;

end

end %--END OF FUNCTION


%--Extra function needed to use parfor (because parfor can have only a
%  single changing input argument).
function [I00,varargout] = model_full_norm_wrapper(li, mp)
    modvar = ModelVariables;
    modvar.sbpIndex = mp.full.indsLambdaMat(li,1);
    modvar.wpsbpIndex = mp.full.indsLambdaMat(li,2);
    modvar.zernIndex = 1;
    modvar.starIndex = 1;
    modvar.whichSource = 'star';
    
    if mp.flagFiber
        [Etemp,Efibertemp] = model_full(mp, modvar,'getNorm');
        I00fiber = abs(Efibertemp).^2;
        varargout{1} = I00fiber;
    else
        Etemp = model_full(mp, modvar,'getNorm');
    end
    I00 = max(max(abs(Etemp).^2));
end
