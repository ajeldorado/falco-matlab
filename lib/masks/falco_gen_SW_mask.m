% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate binary (0-1) software masks for the focal plane. 
% This can be used as a field stop, 
% or for making the scoring and correction regions in the focal plane.
%
%--INPUTS:
% inputs: structure with several fields:
%  -pixresFP: pixels per lambda_c/D
%  -rhoInner: radius of inner FPM amplitude spot (in lambda_c/D)
%  -rhoOuter: radius of outer opaque FPM ring (in lambda_c/D)
%  -angDeg: angular opening (degrees) on the left/right/both sides.
%  -whichSide: which sides to have open. 'left','right', 'top', 'bottom', or 'both'
%  -centering: centering of the coordinates. 'pixel' or 'interpixel'
%  -FOV: minimum desired field of view (in lambda_c/D)
%  -shape: 'square' makes a square. Omitting makes a circle. 
%  -clockAngDeg: Dark hole rotation about the z-axis (deg)
%
%--OUTPUTS:
% maskSW: rectangular, even-sized, binary-valued software mask
% xis: vector of coordinates along the horizontal axis (in lambda_c/D)
% etas: : vector of coordinates along the vertical axis (in lambda_c/D)

function [softwareMask, xis, etas] = falco_gen_SW_mask(inputs)     

% REQUIRED USER INPUTS
pixresFP = inputs.pixresFP; %--pixels per lambda_c/D
rhoInner = inputs.rhoInner; % radius of inner FPM amplitude spot (in lambda_c/D)
rhoOuter = inputs.rhoOuter; % radius of outer opaque FPM ring (in lambda_c/D)
angDeg = inputs.angDeg; %--angular opening (input in degrees) on the left/right/both sides of the dark hole.
whichSide = lower(inputs.whichSide); %--which (sides) of the dark hole have open

angRad = angDeg*(pi/180); %--Convert opening angle to radians

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% OPTIONAL USER INPUTS

% minimum +/- field of view along both axes
if isfield(inputs,'FOV') 
    minFOVxi = inputs.FOV;
    minFOVeta = inputs.FOV;
else
    minFOVxi = rhoOuter; %--Default to the size of the viewable area the field of view is not specified.
    minFOVeta = rhoOuter;
end
%--Overwrite FOV values if specified individually
if isfield(inputs,'xiFOV'); minFOVxi = inputs.xiFOV; end % minimum field of view along horizontal (xi) axis
if isfield(inputs,'etaFOV'); minFOVeta = inputs.etaFOV; end % minimum field of view along vertical (eta) axis
    
if isfield(inputs,'centering'); centering = inputs.centering; else; centering = 'pixel'; end %--Default to pixel centering if it is not specified.

 %--shape of the outer part of the dark hole
if isfield(inputs,'shape'); darkHoleShape = inputs.shape; else; darkHoleShape = 'circle'; end %--Default to a circular outer edge

if isfield(inputs,'clockAngDeg');  clockAngDeg = inputs.clockAngDeg; else; clockAngDeg = 0;  end %--Amount extra to clock the dark hole

%--Lateral offsets of the dark hole
if(isfield(inputs,'xiOffset')); xiOffset = inputs.xiOffset; else; xiOffset = 0; end
if(isfield(inputs,'etaOffset')); etaOffset = inputs.etaOffset; else; etaOffset = 0; end

%--Number of points across each axis. Crop the vertical (eta) axis if angDeg<180 degrees.
if( strcmpi(centering,'interpixel') )
    Nxi =  ceil_even(2*minFOVxi*pixresFP); % Number of points across the full FPM
    Neta = ceil_even(2*minFOVeta*pixresFP);
else
    Nxi =  ceil_even(2*(minFOVxi*pixresFP+1/2)); % Number of points across the full FPM
    Neta = ceil_even(2*(minFOVeta*pixresFP+1/2));
end
%--Overwrite the calculated value if it is specified.
if(isfield(inputs, 'Nxi')); Nxi = inputs.Nxi; end
if(isfield(inputs, 'Neta')); Neta = inputs.Neta; end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--Focal Plane Coordinates
dxi = 1/pixresFP;
deta = dxi;
if( strcmpi(centering,'interpixel')  )
    xis  = (-(Nxi-1)/2: (Nxi-1)/2 )*dxi;
    etas = (-(Neta-1)/2:(Neta-1)/2)*deta;
else %--pixel centering
    xis  = (-Nxi/2: (Nxi/2-1) )*dxi;
    etas = (-Neta/2:(Neta/2-1))*deta;
end
[XIS, ETAS] = meshgrid(xis,etas);
XIS = XIS - xiOffset;
ETAS = ETAS - etaOffset;
[THETAS, RHOS] = cart2pol(XIS, ETAS);

if any(strcmp(whichSide, {'r', 'right', 'lr', 'rl', 'leftright', 'rightleft', 'both'}))
    clockAngRad = 0;
elseif any(strcmp(whichSide, {'l', 'left'}))
    clockAngRad = pi;
elseif any(strcmp(whichSide, {'t', 'top', 'u', 'up', 'tb', 'bt', 'ud', 'du', 'topbottom', 'bottomtop', 'updown', 'downup'}))
    clockAngRad = pi/2;
elseif any(strcmp(whichSide, {'b', 'bottom', 'd', 'down'}))
    clockAngRad = 3/2*pi;
else
    error('falco_gen_SW_mask.m: Invalid value given for inputs.whichSide')
end
clockAngRad = clockAngRad + clockAngDeg*pi/180; %--Add extra clocking specified by inputs.clockAngDeg


% if strcmpi(whichSide,'l') || strcmpi(whichSide,'left') || strcmpi(whichSide,'lr') || strcmpi(whichSide,'lr')
%     clockAngRad = pi;
% elseif strcmpi(whichSide,'r') || strcmpi(whichSide,'right')
%     clockAngRad = 0;
% elseif strcmpi(whichSide,'t') || strcmpi(whichSide,'top')
%     clockAngRad = pi/2;
% elseif strcmpi(whichSide,'b') || strcmpi(whichSide,'bottom')
%     clockAngRad = 3/2*pi;   
% elseif strcmpi(whichSide,'both')
%     clockAngRad = 0;
% else
%     error('falco_gen_SW_mask.m: Unknown value specified for inputs.whichSide')
% end
% clockAngRad = clockAngRad + clockAngDeg*pi/180; %--Add extra clocking specified by inputs.clockAngDeg

%--Generate the Outer Mask
switch lower(darkHoleShape)
    case{'circle', 'annulus'}
        softwareMask0 = (RHOS>=rhoInner & RHOS<=rhoOuter);
    case{'square'}
        softwareMask0 = (RHOS>=rhoInner & abs(XIS)<=rhoOuter & abs(ETAS)<=rhoOuter);
    case{'rect', 'rectangle'}
        softwareMask0 = (RHOS.*cos(THETAS-clockAngRad)>=rhoInner & RHOS.*cos(THETAS-clockAngRad)<=rhoOuter & RHOS.*sin(THETAS-clockAngRad)<=rhoOuter & RHOS.*sin(THETAS-clockAngRad)>=-rhoOuter) |...
                        (RHOS.*cos(THETAS-clockAngRad)<=-rhoInner & RHOS.*cos(THETAS-clockAngRad)>=-rhoOuter & RHOS.*sin(THETAS-clockAngRad)<=rhoOuter & RHOS.*sin(THETAS-clockAngRad)>=-rhoOuter);
    case{'d'}
        softwareMask0 = ((RHOS.*cos(THETAS-clockAngRad)>=rhoInner | RHOS.*cos(THETAS-clockAngRad)<=-rhoInner) & RHOS<=rhoOuter);
    otherwise
        error('falco_gen_SW_mask.m: Invalid value given for inputs.shape')
end

softwareMask = softwareMask0 & abs(angle(exp(1i*(THETAS-clockAngRad))))<=angRad/2;

if any(strcmpi(whichSide, {'both', 'lr', 'rl', 'leftright', 'rightleft', 'tb', 'bt', 'ud', 'du', 'topbottom', 'bottomtop', 'updown', 'downup'}))
    softwareMask2 = softwareMask0 & abs(angle(exp(1i*(THETAS-(clockAngRad+pi)))))<=angRad/2;
    softwareMask = or(softwareMask, softwareMask2);
end

end %--END OF FUNCTION