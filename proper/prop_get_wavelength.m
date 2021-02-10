%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function wlm = prop_get_wavelength(bm)
%        wlm = prop_get_wavelength(bm)
% Return the wavelength of the current beam.
%
% Outputs:
% wlm  = current beam's wavelength (m)
%
% Required inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  wlm  = bm.wl;
end                     % function prop_get_wavelength
