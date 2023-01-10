% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute a radial butterworth filter in 2-D.
%
%--INPUTS
%  nrow: number of rows for the 2-D output array
%  ncol: number of columns for the 2-D output array
%  exponent: exponent of the filter
%  D0: 
%
%--OUTPUTS
%  butter2d: 2-D map containing the butterworth filter
%

function butter2d = falco_butterworth(nrow, ncol, exponent, D0)

    u = 0:nrow-1;
    v = 0:ncol-1;
    idx = find(u>nrow/2);
    idy = find(v>ncol/2);
    u(idx) = u(idx) - nrow;
    v(idy) = v(idy) - ncol;
    [V, U] = meshgrid(v, u);
    D = sqrt(U.^2 + V.^2);
    butter2d = fftshift(1./(1+(D./D0).^(2*exponent)));

end %--END OF FUNCTION
