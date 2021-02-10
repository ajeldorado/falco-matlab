%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wavefront, sampling] = multi_example(lambda_m, n, optval)
%        [wavefront, sampling] = multi_example(lambda_m, n, optval)
% A multi_example system from the Proper Manual
%
% Outputs:
% wavefront  = 2D wavefront array complex E field
% sampling   = sampling distance (m)
%
% Required inputs:
% lambda_m   = wavelength (m)
% n          = number of pixels
%
% Optional inputs:
% optval     = structure of optional values

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
% 2017 Nov 14  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  diam        =    0.048d0   ;  % diameter (m)
  fl_lens     =    0.480d0   ;  % focal length of lens (m)
  n_actuators =   48         ;  % number of DM actuators in each dimension
  pupil_ratio =    0.250d0   ;  % initial beam width/grid width

  wavefront = prop_begin(diam, lambda_m, n, pupil_ratio);
  wavefront = prop_circular_aperture(wavefront, diam / 2.0d0);
  wavefront = prop_define_entrance(wavefront);
  if (nargin > 2) & (optval.use_dm == 1)
    dm_xc = fix(n_actuators / 2.0);     % index of center actuator X
    dm_yc = fix(n_actuators / 2.0);     % index of center actuator Y
    dm_spacing =  1.0d-3;               % DM actuator spacing (m)
    wavefront = prop_dm(wavefront, optval.dm, dm_xc, dm_yc, dm_spacing);
  end
  wavefront = prop_lens(wavefront, fl_lens);

  wavefront = prop_propagate(wavefront, fl_lens);

  [wavefront, sampling] = prop_end(wavefront, 'noabs');

end                     % function multi_example
