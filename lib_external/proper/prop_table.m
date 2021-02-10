%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function prop_table
%        prop_table
% Print the beam dimensions and sampling at each surface in a prescription.
% Intended for internal use by PROPER routines.

% 2005 Feb     jek  created idl routine
% 2015 Jun 22  gmg  Matlab translation
% 2017 Apr 11  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  fprintf(1, '                          Dist from     ');
  fprintf(1, '                Paraxial\n');

  fprintf(1, '                           previous     ');
  fprintf(1, '  Sampling/       beam           Grid   ');
  fprintf(1, '    Lens paraxial     Paraxial\n');

  fprintf(1, 'Surface    Surface         surface      ');
  fprintf(1, '    pixel       diameter         width  ');
  fprintf(1, '    focal length       focal\n');

  fprintf(1, '  type       name            (m)        ');
  fprintf(1, '    (m)            (m)            (m)   ');
  fprintf(1, '         (m)           ratio\n');

  fprintf(1, '--------  ------------  -------------   ');
  fprintf(1, '------------   ------------   ----------');
  fprintf(1, '--  -------------   ------------\n');

  propcommon

  for il = 1 : ActionNum
    if fll(il) ~= 0
    fprintf(1, '   LENS  %12s  %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n', ...
    snml{il}, dzl(il), saml(il), bmdl(il), npa * saml(il), fll(il), efrl(il));

    else
    fprintf(1, 'SURFACE  %12s  %14.6e %14.6e %14.6e %14.6e\n', ...
    snml{il}, dzl(il), saml(il), bmdl(il), npa * saml(il));

    end
  end
end                     % function prop_table
