%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function da = prop_get_sampling_arcsec(bm)
%        da = prop_get_sampling_arcsec(bm)
% Return the current wavefront sampling interval in arcseconds / pixel
% This routine is only valid when the current wavefront is at focus.
%
% Outputs:
% da   = sampling interval (arcseconds)
%
% Required inputs:
% bm   = beam structure

% 2005 Feb     jek  created idl routine
% 2014 May 07  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  fl   = prop_get_fratio(bm) * bm.diam;
  da   = prop_get_sampling(bm) * 360. * 3600. / (2. * pi * fl);
end                     % function prop_get_sampling_arcsec
