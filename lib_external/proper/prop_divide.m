%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_divide(bm, div)
%        bm = prop_divide(bm, div)
% Divides the wavefront array by a value or an array.
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)
% div  = either a 2D array (the same size as the wavefront array) or a scalar
%        The div array is assumed to be centered at pixel(ny/2+1, nx/2+1).

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation
% 2017 Mar 22  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  [my, mx] = size(div);
  if my * mx > 1                % div is a 2D array
    [ny, nx] = size(bm.wf);
    if mx ~= nx | my ~= ny      % div size does not equal wavefront size
      error('Proper:PROP_DIVIDE', ...
      'div array size (%d, %d) not equal to wf array size (%d, %d).\n', ...
      my, mx, ny, nx);
    end
    bm.wf = bm.wf ./ prop_shift_center(div, 'inv');
  else
    bm.wf = bm.wf / div;
  end
end                     % function prop_divide
