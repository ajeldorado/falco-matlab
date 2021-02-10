%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function run_example(wavelength, gridsize)
%        run_example(wavelength, gridsize)
% A function to run example_system
%
% Required inputs:
% wavelength  = wavelength (um)
% gridsize    = number of points

% 2005 Feb     jek  created idl routine
% 2017 Feb 08  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
% 2017 Oct 12  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  prop_init_savestate;

  for it = 0 : 10
    psf  = prop_run('example_system', wavelength, gridsize);

% Let us pretend that we now do something useful with
% this iteration's PSF and then compute another.
  end

  prop_end_savestate;

end                     % function run_example
