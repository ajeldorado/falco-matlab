%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_define_entrance(bm)
%        bm = prop_define_entrance(bm)
% Normalize the wavefront amplitude for a total power of 1.
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)

% 2005 Feb     jek  created idl routine
% 2014 Jun 10  gmg  Matlab translation
% 2017 Mar 22  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

  prop_total_intensity  = sum(sum(abs(bm.wf).^2));
  bm.wf = bm.wf / sqrt(prop_total_intensity);
end                     % function prop_define entrance
