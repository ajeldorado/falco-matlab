%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wavefront, sampling] = example_system(wavelength, gridsize)
%        [wavefront, sampling] = example_system(wavelength, gridsize)
% An example system from the Proper Manual
%
% Outputs:
% wavefront   = 2D wavefront array of intensity
% sampling    = sampling distance (m)
%
% Required inputs:
% wavelength  = wavelength (m)
% gridsize    = number of pixels

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
% 2017 Oct 12  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  diam        =    1.000d0      ;  % diameter (m)
  lens_fl     =   20.000d0      ;  % focal length of each lens (m)
  beam_ratio  =    0.500d0      ;  % initial beam width/grid width

  wavefront   = prop_begin(diam, wavelength, gridsize, beam_ratio);

  if (prop_is_statesaved(wavefront) == 0)
    wavefront   = prop_circular_aperture(wavefront, diam / 2.0d0);
    wavefront   = prop_define_entrance(wavefront);
    wavefront   = prop_lens(wavefront, lens_fl, '1st lens');
    wavefront   = prop_propagate(wavefront, lens_fl, ...
                    'surface_name', 'intermediate focus');
  end

  wavefront   = prop_state(wavefront);

% We are now at the intermediate focus, so pretend that
% we do something to the wavefront here and continue on.

  wavefront   = prop_propagate(wavefront, lens_fl, ...
                  'surface_name', 'second lens');
  wavefront   = prop_lens(wavefront, lens_fl, 'second lens');

  [wavefront, sampling] = prop_end(wavefront);

end                     % function example_system
