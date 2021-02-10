%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_add_wavefront(bm, wfad)
%        bm = prop_add_wavefront(bm, wfad)
% Add a wavefront to the current wavefront array.
% The wavefront array is assumed to be at the same sampling as the
% current wavefront.  Note that this is wavefront, not surface, error.
%
% Outputs:
% bm   = wavefront structure (output)
%
% Required inputs:
% bm   = wavefront structure (input)
% wfad = scalar or 2D image containing the value or wavefront to add.
%        NOTE: All responsibility is on the user to ensure that the
%        two fields have the same sampling and reference phase curvature.
%        NO CHECKING is done by this routine.

% 2011 Mar     jek  created idl routine
% 2016 Apr 05  gmg  Matlab translation
% 2017 Mar 21  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if size(wfad, 1) == 1
    bm.wf = bm.wf + wfad;
  else
    bm.wf = bm.wf + prop_shift_center(wfad, 1);
  end
end                     % function prop_add_wavefront
