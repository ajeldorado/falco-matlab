%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_add_phase(bm, perr)
%        bm = prop_add_phase(bm, perr)
% Add a phase error map or value to the current wavefront array.
% The phase error array is assumed to be at the same sampling as the
% current wavefront.  Note that this is wavefront, not surface, error.
%
% Outputs:
% bm   = wavefront structure (output)
%
% Required inputs:
% bm   = wavefront structure (input)
% perr = scalar or 2D image containing phase error in meters

% 2005 Feb     jek  created idl routine
% 2014 May 08  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed prop_shift_center call to allow odd size arrays
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if size(perr, 1) == 1
    bm.wf = bm.wf .* exp(2.0 * pi * i * perr / bm.wl);
  else
    bm.wf = bm.wf .* exp(2.0 * pi * i * prop_shift_center(perr, 'inv') / bm.wl);
  end
end                     % function prop_add_phase
