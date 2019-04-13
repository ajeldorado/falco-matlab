% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute the next highest even integer of the input.
%
%--INPUTS
% x_in: Scalar value.
%
%--OUTPUTS
%  x_out: even-valued integer value
%
%--VERSION CHANGE HISTORY
% Written on 2017-11-17 by A.J. Riggs.

function x_out = ceil_even(x_in)

    x_out = 2*ceil( 0.5*x_in );

end %--END OF FUNCTION