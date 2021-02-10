%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function dr = prop_get_sampling_radians(bm)
%        dr = prop_get_sampling_radians(bm)
% Return the current wavefront sampling interval in radians / pixel
% This routine is only valid when the current wavefront is at focus.
%
% Outputs:
% dr   = sampling interval (radians)
%
% Required inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  fl   = prop_get_fratio(bm) * bm.diam;
  dr   = prop_get_sampling(bm) / fl;
end                     % function prop_get_sampling_radians
