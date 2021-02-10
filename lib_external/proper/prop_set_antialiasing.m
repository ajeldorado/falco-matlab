%   Copyright 2020 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology


function prop_set_antialiasing( nsub )
%
% This function sets the pixel subsampling factor (nsub by nsub) used to 
% antialias the edges of shapes.  NOTE: Must be called AFTER prop_begin.
%
% Required inputs:
%   nsub: Subsampling factor (must be odd-valued integer)
% 
% 2020 Jan 29  jek  created routine
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

  s = int16(nsub);
  if mod(s,2) == 0
    error(':PROP_SET_ANTIALIASING', 'Subsampling factor must be odd-valued integer.');
  end

  antialias_subsampling = s;

end
