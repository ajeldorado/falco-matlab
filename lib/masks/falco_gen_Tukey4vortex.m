% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% REQUIRED INPUTS: 
% Nwindow =
% RHO = 
% alpha  = 
%
% OUTPUTS:
%  w:     2-D square array of the specified Tukey window
%
% Written by Garreth Ruane.

function w = falco_gen_Tukey4vortex( Nwindow, RHO, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nlut = round(10*Nwindow);
p = linspace(-Nwindow/2,Nwindow/2,Nlut);
lut = tukeywindow(Nlut,alpha);

w = interp1(p,lut,RHO,'linear',0);
end

