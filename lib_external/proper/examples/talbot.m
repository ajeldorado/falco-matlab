%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wavefront, sampling] = talbot(wavelength, gridsize, optval)
%        [wavefront, sampling] = talbot(wavelength, gridsize, optval)
% The Talbot Effect example from the Proper Manual
%
% Outputs:
% wavefront     = 2D wavefront complex array
% sampling      = pixel sampling (m)
%
% Required inputs:
% wavelength    = wavelength (m)
% gridsize      = number of pixels
% optval        = optional values structure
% optval.diam   = diameter (m)
% optval.dist   = propagation distance (m)
% optval.period = period (m)

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
% 2017 Nov 13  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  talbot_length = 2.0d0 * optval.period^2 / wavelength; % Talbot length (m)

  wavefront = prop_begin(optval.diam, wavelength, gridsize);

% Create 1-D grating pattern
  m    =    0.200d0          ;  % pattern amplitude
  x    = ((1:gridsize)-fix(gridsize/2)-1)*prop_get_sampling(wavefront);

% Create 2-D amplitude grating pattern
  [gx, gy] = meshgrid(x, x);
  grating  = 0.5d0 * (1.0d0 + m * cos(2.0d0*pi*gx/optval.period));

  wavefront = prop_multiply(wavefront, grating);
  wavefront = prop_define_entrance(wavefront);
  wavefront = prop_propagate(wavefront, optval.dist, 'to_plane');
  [wavefront, sampling] = prop_end(wavefront, 'noabs');

end                     % function talbot
