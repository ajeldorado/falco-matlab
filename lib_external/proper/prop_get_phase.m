%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function ph = prop_get_phase(bm)
%        ph = prop_get_phase(bm)
% Return the phase of the wavefront array
%
% Outputs:
% ph   = phase of the wavefront array
%
% required inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ph   = prop_shift_center(angle(bm.wf));
end                     % function prop_get_phase
