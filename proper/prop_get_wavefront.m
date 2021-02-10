%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function wf = prop_get_wavefront(bm)
%        wf = prop_get_wavefront(bm)
% Return the complex-valued wavefront array
%
% Outputs:
% wf   = complex-valued wavefront array
%
% Required inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2014 Jul 29  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  wf   = prop_shift_center(bm.wf);
end                     % function prop_get_wavefront
