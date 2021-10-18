%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wavefront, sampling] = psdtest(wavelength, gridsize, usepsdmap)
%        [wavefront, sampling] = psdtest(wavelength, gridsize, usepsdmap)
%
% Outputs:
% wavefront  = 2D wave front array intensity
% sampling   = 2D wave front array sampling distance (m)
%
% Required inputs:
% wavelength = wavelength (m)
% gridsize   = number of pixels
%
% Optional inputs:
% usepsdmap           : if not set, read in WFE map in errormap.fits
%                     : if set, generate and use psd-defined map

% 2005 Feb     jek  created idl routine
% 2017 Apr 21  gmg  Matlab translation
% 2017 Oct 31  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Default parameters
  if nargin < 3
    usepsdmap =    0         ;  % Read in WFE map in errormap.fits
  end

  beam_width_ratio =    0.500d0   ; % beam diameter fraction
  lens_diam        =    0.212d0   ; % telescope diameter (m)
  lens_fl =   24.000d0 * lens_diam; % focal length (m) (f/24 focal ratio)

  wavefront = prop_begin(lens_diam, wavelength, gridsize, beam_width_ratio);

% Create circular entrance aperture
  wavefront = prop_circular_aperture(wavefront, lens_diam / 2.0d0);
  wavefront = prop_define_entrance(wavefront);

% If the variable usepsmap is not set, read in and use the map which
% represents the wavefront error (in this case, it's in nanometers,
% hence we need to multiply it by 1e-9 to convert it to meters) and has
% a sampling of 0.4 mm/pixel.  If usepsdmap is set, then generate and
% use a PSD-defined map.  The maps have an RMS of about 1.0 nm.

  if usepsdmap == 0
    wavefront = prop_errormap(wavefront, 'errormap.fits', 'wavefront', ...
      'sampling', 0.0004d0, 'multiply', 1.000d-9)
  else
    a =    3.290d-23  ; % low spatial frequency power (m^4)
    b =  212.260d0    ; % correlation length parameter (cycles / meter)
    c =    7.800d0    ; % high frequency falloff power law exponent
    wavefront = prop_psd_errormap(wavefront, a, b, c);
  end

  wavefront = prop_lens(wavefront, lens_fl, 'telescope lens');
  dif = prop_get_distancetofocus(wavefront);    % distance to focus (m)
  wavefront = prop_propagate(wavefront, dif, 'snm', 'intermediate focus');

% Multiply field by occulting mask with 4*wavelength/lens_diam HWHM transmission

  [wavefront, mask] = prop_8th_order_mask(wavefront, 4.0d0, 'circular');

  wavefront = prop_propagate(wavefront, lens_fl, 'snm', 'pupil imaging lens');
  wavefront = prop_lens(wavefront, lens_fl, 'pupil imaging lens');
  wavefront = prop_propagate(wavefront, lens_fl, 'snm', 'lyot stop');
  wavefront = prop_circular_aperture(wavefront, 0.53d0, 'norm');  % Lyot stop

  wavefront = prop_propagate(wavefront, lens_fl, 'snm', 'reimaging lens');
  wavefront = prop_lens(wavefront, lens_fl, 'reimaging lens');
  dff = prop_get_distancetofocus(wavefront);    % distance to focus (m)
  wavefront = prop_propagate(wavefront, dff, 'snm', 'final focus');

  [wavefront, sampling] = prop_end(wavefront);

end                     % function psdtest
