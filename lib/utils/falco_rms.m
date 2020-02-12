% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute the RMS of an array.
%  This is needed because Matlab's rms() function is part of the Signal
%  Processing Toolbox.
%
%--INPUTS
%  x: array of numbers
%
%--OUTPUTS
%  RMS: the root-mean-square of the numerical array
%
%--VERSION CHANGE HISTORY
% Written on 2019-02-25 by A.J. Riggs.

function RMS = falco_rms(x)

    RMS = sqrt( mean( abs(x(:)).^2 ) );

end %--END OF FUNCTION