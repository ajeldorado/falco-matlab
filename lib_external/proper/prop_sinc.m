%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function sncu = prop_sinc(xrad)
%        sncu = prop_sinc(xrad)
%
% Outputs:
% sncu = unnormalized sinc function
%
% Required inputs:
% xrad = sncu argument (radians)

% Revision history:
% 2005 Feb     jek  created idl routine
% 2016 Jun 22  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  sncu = ones(size(xrad));
  sncu(xrad ~= 0) = sin(xrad(xrad ~= 0)) ./ xrad(xrad ~= 0);

end                     % function prop_sinc
