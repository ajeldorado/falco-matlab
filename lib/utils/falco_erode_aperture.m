% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ outputAP ] = falco_erode_aperture( inputAP, pixelsInSE )
%falco_erode_aperture Uses imerode to undersize the aperture stop with arbitrary
%shape. 
%   Detailed explanation goes here

% Define the structuring element
se = strel('disk',pixelsInSE); % disk, radius 1

% Erode the rounded aperture
outputAP = imerode(round(inputAP),se);

end