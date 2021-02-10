%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function dist = prop_get_distancetofocus(bm)
%        dist = prop_get_distancetofocus(bm)
% Determine the distance to focus from the current location (m)
%
% Outputs:
% dist = distance to focus (m)
%
% Required inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2014 Jul 29  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dist = bm.w0_pz - bm.pz;
end                     % function prop_get_distancetofocus
