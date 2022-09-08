% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% [sinMode, cosMode] = falco_make_fourier_mode(Nact, NactBeam, xi, eta)
% 
% Compute the DM commands to make a sine wave and a cosine wave at the
% given spatial frequency.
%
% INPUTS
% Nact: Number of actuators across the square DM grid.
% NactBeam: Number of actuators across illuminated beam. Used to rescale
%   the coordinates so that spatial frequencies correspond to the beam
%   diameter instead of the DM width.
% xi: horizontal axis location of the spatial frequency
% eta: vertical axis location of the spatial frequency
%
% OUTPUTS
% sinMode: 2-D array of the sine mode for the given spatial frequency
% cosMode: 2-D array of the cosine mode for the given spatial frequency

function [sinMode, cosMode] = falco_make_fourier_mode(Nact, NactBeam, xi, eta)

    xis = (-(Nact-1)/2:(Nact-1)/2)/NactBeam;  % centered on 0 for even or odd values
    [XIS, ETAS] = meshgrid(xis);
    
    sinMode = sin(2*pi*(xi*XIS + eta*ETAS));
    cosMode = cos(2*pi*(xi*XIS + eta*ETAS));

end %--END OF FUNCTION
