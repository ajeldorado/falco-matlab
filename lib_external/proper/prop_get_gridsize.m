%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function np = prop_get_gridsize()
%        np = prop_get_gridsize()
% Return the dimensions of the current wavefront
%
% Outputs:
% np   = grid size

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon;
  np   = npa;
end                     % function prop_get_gridsize
