%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function prop_init_savestate
%        prop_init_savestate
% Initialize the save state system.  This must be called
% before any calls to prop_state or prop_is_statesaved.
% This routine should be called only once during a session.

% 2005 Feb     jek  created idl routine
% 2016 Apr 06  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

  save_state = 1;
  save_state_lam = 0;

  tmns = double(tic);   % time (ns)
  num  = fix(tmns - fix(tmns * 1.0d-6) * 1.0d+6);
  statefile = ['_' sprintf('%06d', num) '_prop_savestate'];

end                     % function prop_init_savestate
