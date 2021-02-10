%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function arrs = prop_shift_center(arr, inv)
%        arrs = prop_shift_center(arr, inv)
% Shift an array from the origin (1, 1) to the center (default)
% or from the center to the origin ('inv')
%
% Outputs:
% arrs = complex or real array shifted (output)
%
% Required inputs:
% arr  = complex or real array (input)
%
% Optional inputs:
% 'inv'               : use inverse shift

% 2005 Feb     jek  created idl routine
% 2014 May 08  gmg  Matlab translation
% 2014 Sep 02  gmg  Added inverse shift to correctly handle odd size arrays
% 2017 Feb 28  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if nargin == 2 & inv == 'inv'
    arrs = ifftshift(arr);
  else
    arrs =  fftshift(arr);
  end
end                     % function prop_shift_center
