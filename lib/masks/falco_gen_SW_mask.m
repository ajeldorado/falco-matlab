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

%--Read in user inputs
pixresFP = inputs.pixresFP; %--pixels per lambda_c/D
rhoInner = inputs.rhoInner; % radius of inner FPM amplitude spot (in lambda_c/D)
rhoOuter = inputs.rhoOuter; % radius of outer opaque FPM ring (in lambda_c/D)
angDeg = inputs.angDeg; %--angular opening (input in degrees) on the left/right/both sides of the dark hole.
whichSide = inputs.whichSide; %--which (sides) of the dark hole have open

if( isfield(inputs,'FOV') ) % minimum field of view along horizontal (xi) axis
    minFOV = inputs.FOV;
else
    minFOV = rhoOuter; %--Default to the size of the viewable area the field of view is not specified.
end

if( isfield(inputs,'centering') )
    centering = inputs.centering;
else
    centering = 'pixel'; %--Default to pixel centering if it is not specified.
end

if( isfield(inputs,'shape') ) %--shape of the outer part of the dark hole
    darkHoleShape = inputs.shape;
else
    darkHoleShape = 'circle'; %--Default to a circular outer edge
end

if(~isfield(inputs,'clockAngDeg'));  inputs.clockAngDeg = 0;  end %--Amount extra to clock the dark hole


xiCenter = 0;
etaCenter = 0;

if(isfield(inputs,'xi_cen'))
    xiCenter = inputs.xi_cen;
end

if(isfield(inputs,'eta_cen'))
    etaCenter = inputs.eta_cen;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--Convert opening angle to radians
angRad = angDeg*(pi/180);

%--Number of points across each axis. Crop the vertical (eta) axis if angDeg<180 degrees.
if( strcmpi(centering,'interpixel') )
    Nxi =  ceil_even(2*minFOV*pixresFP); % Number of points across the full FPM
    Neta = ceil_even(2*minFOV*pixresFP);
else
    Nxi =  ceil_even(2*(minFOV*pixresFP+1/2)); % Number of points across the full FPM
    Neta = ceil_even(2*(minFOV*pixresFP+1/2));
end

%--Overwrite the calculated value if it is specified.
if(isfield(inputs,'Nxi'))
    Nxi = inputs.Nxi;
end

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
[XIS,ETAS] = meshgrid(xis,etas);
% [THETA,RHO] = cart2pol(XIS,ETAS);
THETAS = atan2((XIS-xiCenter), (ETAS-etaCenter));
RHOS = sqrt((XIS - xiCenter).^2 + (ETAS - etaCenter).^2);

%--Generate the Outer Mask
switch lower(darkHoleShape)
    case{'square'}
        softwareMask0 = (RHOS>=rhoInner & abs(XIS)<=rhoOuter & abs(ETAS)<=rhoOuter);
    otherwise
        softwareMask0 = (RHOS>=rhoInner & RHOS<=rhoOuter);
end

if( strcmpi(whichSide,'L') || strcmpi(whichSide,'left') )
    clockAngRad = 3*pi/2;
elseif( strcmpi(whichSide,'R') || strcmpi(whichSide,'right') )
    clockAngRad = pi/2;
elseif( strcmpi(whichSide,'T') || strcmpi(whichSide,'top') )
    clockAngRad = 0;
elseif( strcmpi(whichSide,'B') || strcmpi(whichSide,'bottom') )
    clockAngRad = pi;   
elseif( strcmpi(whichSide,'both') )
    clockAngRad = 0;
else
    error('falco_gen_SW_mask.m: Unknown value specified for inputs.whichSide')
end

clockAngRad = clockAngRad + inputs.clockAngDeg*pi/180; %--Add extra clocking specified by inputs.clockAngDeg

softwareMask = softwareMask0 & abs(angle(exp(1i*(THETAS-clockAngRad))))<=angRad/2;

if(strcmpi(whichSide,'both'))
    softwareMask2 = softwareMask0 & abs(angle(exp(1i*(THETAS-(clockAngRad+pi)))))<=angRad/2;
    softwareMask = or(softwareMask, softwareMask2);
end

end %--END OF FUNCTION