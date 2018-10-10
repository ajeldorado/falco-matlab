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

function mask = falco_gen_annular_FPM(inputs,varargin)

% Set default values of input parameters
flagRot180deg = false;
%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'rot180'}
            flagRot180deg = true; % For even arrays, beam center is in between pixels.
        otherwise
            error('falco_gen_annular_FPM: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

% %--DEBUGGING ONLY: HARD-CODED INPUTS
% clear all
% addpath ~/Repos/FALCO/lib/PROPER
% inputs.pixresFPM = 6; %--pixels per lambda_c/D
% inputs.rhoInner = 2.8; % radius of inner FPM amplitude spot (in lambda_c/D)
% inputs.rhoOuter = 10.1; % radius of outer opaque FPM ring (in lambda_c/D)
% inputs.FPMampFac = 0.4; % amplitude transmission of inner FPM spot
% inputs.centering = 'interpixel';
% flagRot180deg = false;


pixresFPM = inputs.pixresFPM; %--pixels per lambda_c/D
rhoInner = inputs.rhoInner; % radius of inner FPM amplitude spot (in lambda_c/D)
rhoOuter = inputs.rhoOuter; % radius of outer opaque FPM ring (in lambda_c/D)
FPMampFac = inputs.FPMampFac; % amplitude transmission of inner FPM spot
centering = inputs.centering; % Centering of array: 'pixel' or 'interpixel'

dx = 1/pixresFPM; %--lambda_c/D per pixel.

% % Centering of array: 'pixel' or 'interpixel'
% if(isfield(inputs,'centering'))
%     centering = inputs.centering;
% else %--Default to pixel centering
%     centering = 'pixel';
% end

if(rhoOuter==inf)
    if(strcmpi(centering,'interpixel'))
        Narray = ceil_even((2*rhoInner/dx)); % number of points across the inner diameter of the FPM.
    else
        Narray = ceil_even(2*(rhoInner/dx+1/2)); % number of points across the inner diameter of the FPM. Another half pixel added for pixel-centered masks.
    end
    
    Dmask = 2*pixresFPM*rhoInner; %--Diameter of the mask

else
    if(strcmpi(centering,'interpixel'))
        Narray = ceil_even(2*rhoOuter/dx); % number of points across the outer diameter of the FPM. 
    else
        Narray = ceil_even(2*(rhoOuter/dx+1/2)); % number of points across the outer diameter of the FPM. Another half pixel added for pixel-centered masks.
    end
    
    Dmask = 2*pixresFPM*rhoOuter; %--Diameter of the mask

end

xshift = 0;%inputs.xshift; % translation in x of FPM (in lambda_c/D)
yshift = 0;%inputs.yshift; % translation in y of FPM (in lambda_c/D)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--GENERAL SPECS
Darray = Narray*dx; %--width of array in lambda_c/D
wl_dummy   = 1e-6;              % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 for pixel-centered FPM, or -diam/Narray for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
        if(flagRot180deg)
            cshift = -dx;
        end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--INITIALIZE PROPER. Note that:  bm.dx = diam / bdf / np;
bdf = Dmask/Darray; %--beam diameter factor in output array
bm = prop_begin(Dmask, wl_dummy, Narray,'beam_diam_fraction',bdf);
% figure(1); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;

if(rhoOuter ~= inf)
    %--Outer opaque ring of FPM
    ra_OD = (rhoOuter);
    cx_OD = 0 + cshift + xshift;
    cy_OD = 0 + cshift + yshift;
    bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);%, cx, cy, norm);
    % figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;
    % figure(3); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;
end

%--Inner spot of FPM (Amplitude transmission can be nonzero)
ra_ID = (rhoInner); 
cx_ID = 0 + cshift + xshift;
cy_ID = 0 + cshift + yshift;
innerSpot = prop_ellipse(bm, ra_ID,ra_ID,'cx',cx_ID,'cy',cy_ID,'DARK')*(1-FPMampFac) + FPMampFac;
% figure(4); imagesc((abs(innerSpot))); axis xy equal tight; colorbar;

mask = ifftshift(abs(bm.wf)); %--undo PROPER's fftshift
mask = mask.*innerSpot; %--Include the inner FPM spot
% figure(19); imagesc(mask); axis xy equal tight; colorbar;


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










