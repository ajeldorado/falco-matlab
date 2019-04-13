% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function y = sinc(x)

x0 = x;
x(x0==0)= 1;                           
y = sin(pi*x)./(pi*x);                                                     
y(x0==0) = 1;