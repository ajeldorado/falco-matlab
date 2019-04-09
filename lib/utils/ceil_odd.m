% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute the next highest odd integer of the input.
%
%--INPUTS
% x_in: Scalar value.
%
%--OUTPUTS
%  x_out: even-valued integer value
%
%--VERSION CHANGE HISTORY
% Written on 2018-03-03 by A.J. Riggs.

function x_out = ceil_odd(x_in)

    x_out = ceil(x_in);
    if(mod(x_out,2)==0)
        x_out = x_out + 1;
    end

end %--END OF FUNCTION