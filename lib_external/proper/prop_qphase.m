%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_qphase(bm, cu)
%        bm = prop_qphase(bm, cu)
% This routine applies a quadratic, radially-dependent phase factor
% caused by a wavefront curvature to the current wavefront.
% It is used by prop_stw and prop_wts.
% Intended only for internal use by prop_* routines.
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)
% cu   = phase curvature radius (m)

% 2005 Feb     jek  created idl routine
% 2013 Oct     jek  speed up phase calculation
% 2014 May 12  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed fftshift to ifftshift to allow for odd size arrays
% 2016 Sep 12  gmg  Changed prop_radius2 to prop_radius
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if cu == 0.0
    return
  end

  bm.wf = bm.wf .* exp(i * pi * ifftshift((prop_radius(bm)).^2) / cu / bm.wl);
end                     % function prop_qphase
