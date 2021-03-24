%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function prop_end_savestate
%        prop_end_savestate
% Terminate the current save state system.
% This deletes the files created by prop_state/prop_savestate.
% This routine should be called only once during a session.

% 2005 Feb     jek  created idl routine
% 2016 Apr 06  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

  save_state = 0;
  save_state_lam = 0;

  system(['rm *' statefile]);
end                     % function prop_end_savestate
