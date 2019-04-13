% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = propcustom_2FT(Ein,varargin)
%
%--This function forward propagates 2 Fourier Transforms ahead simply by
%rotating the matrix 180 degrees and maintaining the correct centering of
%the beam on the array.
%
%
%--INPUTS
% Ein = electric field at at input (focal) plane
% 'INTERPIXEL' = optional flag to set the centering of the array between pixels
%
% - Modified on 2019-04-05 by A.J. Riggs to remove the 1/1i term from each FT.
%
%--------------------------------------------------------------------------

function Eout = propcustom_2FT(Ein,varargin)

% Set default values of input parameters
centering = 'pixel';

%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'interpixel'}
            centering = 'interpixel'; % For even arrays, beam center is in between pixels.
        case {'pixel'}
            centering = 'pixel'; % For even arrays, beam center is in between pixels.    
        otherwise
            error('propcustom_2FT: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

Eout = rot90(Ein,2);  %--Forward propagate two Fourier transforms by rotating 180 degrees.

switch centering
    case{'pixel'}
        Eout = circshift(Eout,[1 1]);   %--To undo center offset when beam is pixel centered and rotating by 180 degrees.
    case{'interpixel'}
        %--No change
    otherwise
        error('propcustom_2FT: centering not specified correctly')
end

end %--END OF FUNCTION