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

function mask = falco_gen_annular_FPM_GR(inputs,varargin)

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

    [X,Y] = meshgrid(-Narray/2:Narray/2-1);
    [~,RHO] = cart2pol(X,Y);
    mask = FPMampFac*exp(-(RHO/(rhoInner*pixresFPM)).^500);
    


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









