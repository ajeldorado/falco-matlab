%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function amp = prop_get_amplitude(bm)
%        amp = prop_get_amplitude(bm)
% Return the amplitude of the wavefront array
%
% Outputs:
% amp  = amplitude of the wavefront array
%
% Require inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation
% 2017 Mar 22  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  amp  = prop_shift_center(abs(bm.wf));
end                     % function prop_get_amplitude
