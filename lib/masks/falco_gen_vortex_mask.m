% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% REQUIRED INPUTS: 
% charge  = charge of the vortex
% N   = number of points across the array 
%
% OUTPUTS:
%  V:     2-D square array of the vortex E-field
%
% Created in 2018 by Garreth Ruane.

function V = falco_gen_vortex_mask( charge, N )
%   Detailed explanation goes here

[X,Y] = meshgrid(-N/2:N/2-1);
V = exp(1i*charge*atan2(Y,X));

end