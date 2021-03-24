%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function pz = prop_get_z(bm)
%        pz = prop_get_z(bm)
%
% Outputs:
% pz   = distance from the initialization of the wavefront to the
%        current surface (m)
%
% Required inputs:
% bm   = beam structure

% 2013 Nov     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  pz   = bm.pz;
end                     % function prop_get_z
