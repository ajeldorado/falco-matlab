% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ LS ] = falco_gen_SCDA_LS( LSin, LSout, apRad, EPsupport, LSspiders, Narray, centering  )
%falco_gen_SCDA_LS Summary of this function goes here
%   Detailed explanation goes here

% Defines the coordinate systems
if( strcmp(centering,'interpixel') || strcmp(centering,'even') )
    [X,Y] = meshgrid(-Narray/2+0.5:Narray/2-0.5); % Grids with Cartesian (x,y) coordinates 
else
    [X,Y] = meshgrid(-Narray/2:Narray/2-1); % Grids with Cartesian (x,y) coordinates 
end
[~,RHO] = cart2pol(X,Y);  % Grids with polar (rho,theta) coordinates 

if(LSin > 0)
    LS = exp(-(RHO/(LSout*apRad)).^1000)-exp(-(RHO/(LSin*apRad)).^1000);
else
    LS = exp(-(RHO/(LSout*apRad)).^1000);
end
if(LSspiders > 0)
    LS = ErodeAperture(round(EPsupport),LSspiders).*LS;
end
end

