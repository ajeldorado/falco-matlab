% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = propcustom_relay(Ein,Nrelay,varargin)
%
%--This function relays the input electric field to the next conjugate plane. Each relay 
% causes a 180-degree rotation. The centering has to be 
% corrected after rotation if the beam is pixel-centered on the array.
%
%
%--INPUTS
% Ein = electric field at at input (focal) plane
% Nrelay = number of re-imaging relays to perform
%
%--OPTIONAL INPUTS
% 'INTERPIXEL' or 'PIXEL' = optional flag to set centering of beam in array
%
%--REVISION HISTORY
% - Modified on 2019-04-05 by A.J. Riggs to remove the 1/1i term from each FT.
% - Modified on 2019-02-27 by A.J. Riggs from propcustom_2FT.m to be for an arbitrary 
%   number of re-imaging relays, not just 1.
%--------------------------------------------------------------------------

function Eout = propcustom_relay(Ein,Nrelay,varargin)

    if(Nrelay-round(Nrelay)~=0 || Nrelay<0)
        error('propcustom_relay: The variable Nrelay must be a non-negative integer')
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

    Eout = rot90(Ein,2*Nrelay);  %--Forward propagate two Fourier transforms by rotating 180 degrees. The (1/1j)^2 is from the coefficients of the two FTs.

    if(strcmpi(centering,'pixel') && mod(Nrelay,2)==1)
        Eout = circshift(Eout,[1 1]);   %--To undo center offset when beam is pixel centered and rotating by 180 degrees.
    end

end %--END OF FUNCTION