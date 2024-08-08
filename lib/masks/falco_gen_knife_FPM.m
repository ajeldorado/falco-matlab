% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate an annular FPM in Matlab using PROPER.
%  -Outside the outer ring is opaque.
%   -If rhoOuter = infinity, then the outer ring is omitted and the mask is
%   cropped down to the size of the inner spot.
%  -The inner spot has a specifyable amplitude value.
%  -The output array is the smallest size that fully contains the mask.
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Modified by A.J. Riggs on 2017-10-25 from a file to generate a pupil 
%   to a file that generates an FPM. 
% Modified by A.J. Riggs on 2017-02-09 to be a function with specifyable
%   amounts of magnification, clocking, and translation. 
%   Now keeping the pupil edges non-binary.
% This script was first written by A.J. Riggs on 2017-01-11.
%
%
%--INPUTS:
% pixresFPM: resolution in pixels per lambda_c/D
% rhoInner:   radius of inner FPM amplitude spot (in lambda_c/D)
% rhoOuter:   radius of outer opaque FPM ring (in lambda_c/D). Set to
%             infinity for an occulting-spot only FPM
% amp0:    amplitude transmission of inner FPM spot
%
%--OUTPUT:
% mask:    cropped-down, 2-D FPM representation. amplitude only 

function mask = falco_gen_knife_FPM(inputs)

%--x-offset of mask in lambda0/D
if(isfield(inputs, 'xOffset'))
    xOffset = inputs.xOffset;
else
    xOffset = 0;
end

%--y-offset of mask in lambda0/D
if(isfield(inputs, 'yOffset'))
    yOffset = inputs.yOffset;
else
    yOffset = 0;
end

if(isfield(inputs, 'FPMampFac')); FPMampFac = inputs.FPMampFac; else; FPMampFac = 0; end % amplitude transmission of inner FPM spot
if(isfield(inputs,'centering')); centering = inputs.centering; else; centering = "pixel"; end

pixresFPM = inputs.pixresFPM; %--pixels per lambda_c/D
rhoInner = inputs.rhoInner; % radius of inner FPM amplitude spot (in lambda_c/D)
rhoOuter = inputs.rhoOuter; % radius of outer opaque FPM ring (in lambda_c/D)

dx = 1/pixresFPM; %--lambda_c/D per pixel.
maxAbsOffset = max([abs(xOffset), abs(yOffset)]);
if(isinf(rhoOuter))
    if(strcmpi(centering, 'interpixel'))
        Narray = ceil_even(2*rhoInner/dx + 2*pixresFPM*maxAbsOffset); % number of points across the inner diameter of the FPM.
    elseif(strcmpi(centering, 'pixel'))
        Narray = ceil_even(2*rhoInner/dx + 2*pixresFPM*maxAbsOffset + 1); % number of points across the inner diameter of the FPM. Another half pixel added for pixel-centered masks.
    end
    
    Dmask = 2*pixresFPM*rhoInner; %--Diameter of the mask

else
    if(strcmpi(centering, 'interpixel'))
        Narray = ceil_even(2*rhoOuter/dx + 2*pixresFPM*maxAbsOffset); % number of points across the outer diameter of the FPM. 
    elseif(strcmpi(centering, 'pixel'))
        Narray = ceil_even(2*rhoOuter/dx + 2*pixresFPM*maxAbsOffset + 1); % number of points across the outer diameter of the FPM. Another half pixel added for pixel-centered masks.
    end
    
    Dmask = 2*pixresFPM*rhoOuter; %--Diameter of the mask

end

if(isfield(inputs, 'Narray'))
    Narray = inputs.Narray;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--GENERAL SPECS
Darray = Narray*dx; %--width of array in lambda_c/D
wl_dummy   = 1e-6;              % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 for pixel-centered FPM, or -diam/Narray for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
    otherwise
        error('Invalid value %s for centering',centering)
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--INITIALIZE PROPER. Note that:  bm.dx = diam / bdf / np;
bdf = Dmask/Darray; %--beam diameter factor in output array
bm = prop_begin(Dmask, wl_dummy, Narray,'beam_diam_fraction',bdf);

if(isinf(rhoOuter) == false)
    %--Outer opaque ring of FPM
    ra_OD = (rhoOuter);
    cx_OD = 0 + cshift + xOffset;
    cy_OD = 0 + cshift + yOffset;
    bm = prop_circular_aperture(bm, ra_OD, 'cx', cx_OD, 'cy', cy_OD);%, cx, cy, norm);
end

%--Inner spot of FPM (Amplitude transmission can be nonzero)
ra_ID = (rhoInner); 
cx_ID = 0 + cshift + xOffset;
cy_ID = 0 + cshift + yOffset;
innerSpot = prop_rectangle(bm, ra_ID*2, ra_ID*2, 'cx', cx_ID, 'cy', cy_ID, 'DARK')*(1-FPMampFac) + FPMampFac;

mask = ifftshift(abs(bm.wf)); %--undo PROPER's fftshift
mask = mask.*innerSpot; %--Include the inner FPM spot

end %--END OF FUNCTION

% %--DEBUGGING: Visually verify that mask is centered correctly
% figure(11); imagesc(mask); axis xy equal tight; colorbar; drawnow;
% switch centering 
%     case {'pixel'}
%         maskTemp = mask(2:end,2:end);
%     otherwise
%         maskTemp = mask;
% end
% figure(12); imagesc(maskTemp-rot90(maskTemp,2)); axis xy equal tight; colorbar; 
% title('Centering Check','Fontsize',20); set(gca,'Fontsize',20);
% drawnow;
% 
% sum(sum(mask))