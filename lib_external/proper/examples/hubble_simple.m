%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wavefront, sampling] = hubble_simple(wavelength, grid_n, delta_sec)
%        [wavefront, sampling] = hubble_simple(wavelength, grid_n, delta_sec)
%
% Outputs:
% wavefront   = wavefront intensity 2D array
% sampling    = spacing between wavefront array points (m)
%
% Required inputs:
% wavelength  = wavelength (m)
% grid_n      = number of points
%
% Optional inputs:
% delta_sec   = delta primary to secondary separation (m)

% 2016 Oct 20  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
% 2017 Oct 17  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if nargin < 3
    delta_sec = 0.0d0   ;  % delta primary to secondary spacing (m)
  end

% Hubble telescope optical prescription
  diam           =     2.4d0         ;  % diameter of telescope (m)
  fl_pri         =     5.52085d0     ;  % focal length mirror 1 (m)
  d_pri_sec      =     4.907028205d0 ;  % mirror 1 to mirror 2  (m)
  fl_sec         =    -0.6790325d0   ;  % focal length mirror 2 (m)
  d_sec_to_focus =     6.3919974d0   ;  % mirror 2 to focus     (m)
  or2            =     0.396d0       ;  % mirror 2 obscuration  (m)
  pr             =     0.078d0       ;  % mirror 1 pad radius   (m)
  p1x            =    -0.9066d0      ;  % mirror 1 pad 1 X      (m)
  p1y            =    -0.5538d0      ;  % mirror 1 pad 1 Y      (m)
  p2x            =     0.0d0         ;  % mirror 1 pad 2 X      (m)
  p2y            =     1.0705d0      ;  % mirror 1 pad 2 Y      (m)
  p3x            =     0.9127d0      ;  % mirror 1 pad 3 X      (m)
  p3y            =    -0.5477d0      ;  % mirror 1 pad 3 Y      (m)
  vl             =     2.5d0         ;  % mirror 2 vane length  (m)
  vw             =     0.0264d0      ;  % mirror 2 vane width   (m)
  beam_ratio     =     0.500d0       ;  % beam ratio

  wavefront = prop_begin(diam, wavelength, grid_n, beam_ratio);

% primary mirror aperture
  wavefront = prop_circular_aperture(wavefront, diam / 2.0d0);

% secondary obscuration
  wavefront = prop_circular_obscuration(wavefront, or2);

% secondary vane (vertical)
  wavefront = prop_rectangular_obscuration(wavefront, vw, vl);
% secondary vane (horizontal)
  wavefront = prop_rectangular_obscuration(wavefront, vl, vw);

% primary mirror pad 1
  wavefront = prop_circular_obscuration(wavefront, pr, 'xc', p1x, 'yc', p1y);
% primary mirror pad 2
  wavefront = prop_circular_obscuration(wavefront, pr, 'xc', p2x, 'yc', p2y);
% primary mirror pad 3
  wavefront = prop_circular_obscuration(wavefront, pr, 'xc', p3x, 'yc', p3y);

  wavefront = prop_define_entrance(wavefront);

  wavefront = prop_lens(wavefront, fl_pri, 'primary');    % primary mirror

  wavefront = prop_propagate(wavefront, d_pri_sec + delta_sec, ...
                  'surface_name', 'secondary');

  wavefront = prop_lens(wavefront, fl_sec, 'secondary');  % secondary mirror

  wavefront = prop_propagate(wavefront, d_sec_to_focus + delta_sec, ...
                  'surface_name', 'HST focus');

  [wavefront, sampling] = prop_end(wavefront);

end                     % function hubble_simple
