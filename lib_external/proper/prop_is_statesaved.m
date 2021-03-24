%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function stat = prop_is_statesaved(bm)
%        stat = prop_is_statesaved(bm)
% Determine if a previously saved state exists for the current wavelength.
%
% Outputs:
% stat = status = 1 if a state exists, 0 otherwise
%
% Rquired inputs:
% bm   = current beam structure

% 2005 Feb     jek  created idl routine
% 2016 Apr 06  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

  stat = 0;
  if numel(save_state_lam) ~= 0
    for il = 1 : size(save_state_lam, 1)
% Does a state exist for the current wavelength?
      if save_state_lam(il) == bm.wl
        stat = 1;
      end
    end
  end
end                     % function prop_is_statesaved
