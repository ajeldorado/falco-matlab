% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = propcustom_relay(Ein,Nrelay,varargin)
%
%--This function relays the input electric field to the next conjugate plane.
% Each relay causes a 180-degree rotation. The centering has to be 
% corrected after rotation if the beam is pixel-centered on the array.
%
%--INPUTS
% Ein = electric field at at input (focal) plane
% Nrelay = number of re-imaging relays to perform
%
%--OPTIONAL INPUTS
% 'INTERPIXEL' or 'PIXEL' = optional flag to set centering of beam in array
%--------------------------------------------------------------------------

function Eout = propcustom_relay(Ein, Nrelay, varargin)

    % Error checks on inputs
    if length(size(Ein)) ~= 2
        error('Ein must be a 2-D array.')
    elseif size(Ein, 1) ~= size(Ein, 2)
        error('Ein must be a square array.')
    end
    
    if(Nrelay - round(Nrelay) ~= 0)
        error('propcustom_relay: The variable Nrelay must be an integer')
    end
    
    

    % Set default values of input parameters
    centering = 'pixel';

    %--Look for Optional Keywords
    icav = 0; % index in cell array varargin
    while icav < size(varargin, 2)
        icav = icav + 1;
        switch lower(varargin{icav})
            case {'interpixel'}
                centering = 'interpixel'; % For even arrays, beam center is in between pixels.
            case {'pixel'}
                centering = 'pixel'; % For even arrays, beam center is in between pixels.    
            otherwise
                error('propcustom_relay: Unknown keyword: %s\n', ...
                varargin{icav});
        end
    end

    Eout = rot90(Ein, 2*Nrelay);  %--Forward propagate two Fourier transforms by rotating 180 degrees.

    if strcmpi(centering, 'pixel') && (mod(Nrelay, 2) == 1) && (mod(size(Ein, 1), 2) == 0)
        Eout = circshift(Eout, [1 1]);   %--To undo center offset when beam is pixel centered in even-sized array 
    end

end %--END OF FUNCTION