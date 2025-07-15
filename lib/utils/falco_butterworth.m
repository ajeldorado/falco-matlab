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
%  D0:center 
%
%--OUTPUTS
%  butter2d: 2-D map containing the butterworth filter
%

function butter2d = falco_butterworth(nrow, ncol, ctr_x, ctr_y, exponentx,exponenty, X0, Y0, geometry)

    u = 0:nrow-1;
    v = 0:ncol-1;
    idx = find(u>nrow/2);
    idy = find(v>ncol/2);
    u(idx) = u(idx) - nrow ;
    u = u- ctr_x;
    v(idy) = v(idy) - ncol ;
    v = v - ctr_y;
    [U, V] = meshgrid(u, v);
    D = sqrt(U.^2 + V.^2);
    if strcmp(geometry,'squared')
        butter2d = fftshift(1./(1+(abs(U)./X0).^(2*exponentx))*1./(1+(abs(V)./Y0).^(2*exponenty)));
    else
    %butter2d = butter2d .* fftshift(1./(1+(V./Y0).^(2*exponent)));
    %butter2d = abs(butter2d);
        D0 = (X0+Y0)/2;
        exponent = (exponentx+exponenty)/2;
        butter2d = fftshift(1./(1+(D./D0).^(2*exponent)));
    end

end %--END OF FUNCTION
