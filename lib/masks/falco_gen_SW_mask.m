% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate binary (0-1) software masks for the focal plane. 
% This can be used as a field stop, 
% or for making the scoring and correction regions in the focal plane.
%
% Created on 2018-03-07 by A.J. Riggs
%
%--INPUTS:
% inputs: structure with several fields:
%  -pixresFP: pixels per lambda_c/D
%  -rhoInner: radius of inner FPM amplitude spot (in lambda_c/D)
%  -rhoOuter: radius of outer opaque FPM ring (in lambda_c/D)
%  -angDeg: angular opening (degrees) on the left/right/both sides.
%  -whichSide: which sides to have open. 'left','right', or 'both'
%  -centering: centering of the coordinates. 'pixel' or 'interpixel'
%  -FOV: minimum desired field of view (in lambda_c/D)
%
%--OUTPUTS:
% maskSW: rectangular, even-sized, binary-valued software mask
% xis: vector of coordinates along the horizontal axis (in lambda_c/D)
% etas: : vector of coordinates along the vertical axis (in lambda_c/D)

function [maskSW,xis,etas] = falco_gen_SW_mask(inputs)     

%--Read in user inputs
pixresFP = inputs.pixresFP; %--pixels per lambda_c/D
rhoInner = inputs.rhoInner; % radius of inner FPM amplitude spot (in lambda_c/D)
rhoOuter = inputs.rhoOuter; % radius of outer opaque FPM ring (in lambda_c/D)
angDeg = inputs.angDeg; %--angular opening (input in degrees) on the left/right/both sides of the dark hole.
whichSide = inputs.whichSide; %--which (sides) of the dark hole have open

if( isfield(inputs,'FOV') ) % minimum field of view along horizontal (xi) axis
    FOVmin = inputs.FOV;
else
    FOVmin = rhoOuter; %--Default to the size of the viewable area the field of view is not specified.
end

if( isfield(inputs,'centering') )
    centering = inputs.centering;
else
    centering = 'pixel'; %--Default to pixel centering if it is not specified.
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--Convert opening angle to radians
angRad = angDeg*(pi/180);

%--Number of points across each axis. Crop the vertical (eta) axis if angDeg<180 degrees.
if( strcmpi(centering,'interpixel') )
    Nxi =  ceil_even(2*FOVmin*pixresFP); % Number of points across the full FPM
    Neta = ceil_even(2*sin(angRad/2)*FOVmin*pixresFP);
else
    Nxi =  ceil_even(2*(FOVmin*pixresFP+1/2)); % Number of points across the full FPM
    Neta = ceil_even(2*(sin(angRad/2)*FOVmin*pixresFP+1/2));
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
RHOS = sqrt(XIS.^2 + ETAS.^2);
TAN = atan(ETAS./XIS);

%--Generate the Software Mask
maskSW = (RHOS>=rhoInner & RHOS<=rhoOuter & TAN<=angRad/2 & TAN>=-angRad/2);

%--Determine if it is one-sided or not
if( strcmpi(whichSide,'L') || strcmpi(whichSide,'left') )
    maskSW(XIS>=0) = 0;
elseif( strcmpi(whichSide,'R') || strcmpi(whichSide,'right') )
    maskSW(XIS<=0) = 0;
elseif( strcmpi(whichSide,'T') || strcmpi(whichSide,'top') )
    maskSW(ETAS<=0) = 0;
elseif( strcmpi(whichSide,'B') || strcmpi(whichSide,'bottom') )
    maskSW(ETAS>=0) = 0;    
end



end %--END OF FUNCTION
