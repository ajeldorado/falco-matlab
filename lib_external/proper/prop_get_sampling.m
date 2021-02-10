%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function dx = prop_get_sampling(bm)
%        dx = prop_get_sampling(bm)
% Return the current wavefront sampling interval in meters / pixel
%
% Outputs:
% dx   = sampling interval (m)
%
% Required inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2014 May 07  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dx   = bm.dx;
end                     % function prop_get_sampling
