%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function prop_print_zernikes(numz)
%        prop_print_zernikes(numz)
% Print to the screen the first numz of Noll-ordered Zernike polynomials
% for an unobscured circular aperture.
%
% Required inputs:
% numz = number of Zernike polynomials to print (1 to numz)

% Revision history:
% 2005 Feb     jek
% 2016 May 02  gmg  Translated to Matlab
% 2017 Feb 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  list = prop_noll_zernikes(numz);
  for iz = 1 : numz
    fprintf(1, '%8d  =  %s\n', iz, char(list(iz)));
  end
end                     % prop_print_zernikes
