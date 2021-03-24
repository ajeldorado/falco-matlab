%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function dx = prop_get_nyquistsampling(bm, wl)
%        dx = prop_get_nyquistsampling(bm, wl)
% Return the Nyquist sampling interval for the current wavefront.
% This routine determines the Nyquist sampling interval for
% the current beam, which is the focal_ratio * wavelength / 2.
%
% Outputs:
% dx   = Nyquist sampling interval (m)
%
% Required inputs:
% bm   = beam structure
%
% Optional inputs;
% wl   = wavelength (m)
%        By default, the current wavefront wavelength is used.
%        This parameter can be used when you want to know the Nyquist
%        sampling for a wavelength other than for the current wavefront.

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation
% 2017 Mar 16  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if nargin > 1
    dx   = bm.fr *    wl / 2.0;
  else
    dx   = bm.fr * bm.wl / 2.0;
  end
end                     % function prop_get_nyquistsampling
