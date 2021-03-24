%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wavefront, sampling] = microscope(wavelength, gridsize, focus_offset)
%        [wavefront, sampling] = microscope(wavelength, gridsize, focus_offset)
% A Microscope example from the Proper Manual
%
% Outputs:
% wavefront  = 2D wavefront array of intensity
% sampling   = pixel sampling (m)
%
% Required inputs:
% wavelength = wavelength (m)
% gridsize   = number of pixels
%
% Optional inputs:
% focus_offset = focus_offset (m)

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
% 2017 Oct 10  gmg  Fixed bug in objective diameter
% 2017 Nov 13  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  d_objective   =    0.005d0 ;  % objective diameter (m)
  fl_objective  =    0.010d0 ;  % focal length objective (m)
  fl_eyepiece   =    0.020d0 ;  % focal length eyepiece (m)
  fl_eye        =    0.022d0 ;  % focal length human eye (m)
  beam_ratio    =    0.400d0 ;  % initial beam width/grid width
  d1            =    0.160d0 ;  % standard tube length (m)
  d_intermediate_image = fl_objective + d1  ;  % intermediate image distance (m)

% Compute in-focus distance of object from objective
  d_object  = 1.0d0 / (1.0d0 / fl_objective - 1.0d0 / d_intermediate_image);

  wavefront = prop_begin(d_objective, wavelength, gridsize, beam_ratio);

  wavefront = prop_circular_aperture(wavefront, d_objective / 2.0d0);
  wavefront = prop_define_entrance(wavefront);

% Simulate the diverging wavefront emitted from a point source placed
% "d_object" in front of the objective by using a negative lens (focal
% length = -d_object) placed at the location of the objective.

  if (nargin < 3)
    focus_offset = 0.0d0            ;  % focus offset (m)
  end

  wavefront = prop_lens(wavefront, -(d_object + focus_offset));
  wavefront = prop_lens(wavefront, fl_objective, 'objective');

  wavefront = prop_propagate(wavefront, d_intermediate_image, ...
              'surface_name', 'intermediate image');
  wavefront = prop_propagate(wavefront, fl_eyepiece, ...
              'surface_name', 'eyepiece');
  wavefront = prop_lens(wavefront, fl_eyepiece, 'eyepiece');

  exit_pupil_distance = fl_eyepiece / ...
               (1.0d0 - fl_eyepiece / (d_intermediate_image + fl_eyepiece));
  wavefront = prop_propagate(wavefront, exit_pupil_distance, ...
              'surface_name', 'exit pupil/eye');

  wavefront = prop_lens(wavefront, fl_eye, 'eye');
  wavefront = prop_propagate(wavefront, fl_eye, 'surface_name', 'retina', 'to_plane');

% Extract 2D array of intensity and pixel sampling from beam structure
  [wavefront, sampling] = prop_end(wavefront);

end                     % function microscope
