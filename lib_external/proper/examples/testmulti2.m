%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


% Demo script testmulti2.m

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
% 2017 Nov 14  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  gridsize  = 1024           ;  % number of pixels
  npatterns =    3           ;  % number of patterns
  lambda    =    0.600d0     ;  % wavelength (um)

% Create different Deformable Mirror ripple patterns (50 nm amplitude)

  dm     = zeros(48, 48);
  [x, y] = meshgrid(2.0d0 * pi * [0 : 47] / 47.0d0, ones(1, 48));
  use_dm = 1         ;  % multi_example flag: use deformable mirror
  optval = struct('use_dm', use_dm, 'dm', dm);

  for i = 1 : npatterns
    optval(i).dm = 5.0d-8 * cos(4.0d0 * x * i);
    optval(i).use_dm = 1;
  end

% Generate monochromatic fields in parallel
  [fields,sampling]=prop_run_multi('multi_example',lambda,gridsize,'passvalue',optval);
