%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function fr = prop_get_fratio(bm)
%        fr = prop_get_fratio(bm)
% Return the current wavefront focal ratio 
% This routine computes the current beam's focal ratio by dividing
% the current distance to focus by the current beam diameter.
%
% Outputs:
% fr   = focal ratio
%
% Required inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2014 May 07  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  fr   = bm.fr;
end                     % function prop_get_fratio
