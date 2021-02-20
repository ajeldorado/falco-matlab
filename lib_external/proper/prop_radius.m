%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function radi = prop_radius(bm, norm)
%        radi = prop_radius(bm, norm)
% Returns a 2D array in which the value of each element corresponds to
% the distance of that element from the center of the current wavefront.
% By default, the distance is in meters, unless the norm switch is set,
% in which case it is normalized to the current radius of the beam.
% The center of the wavefront is set to be at the center of the array.
%
% Outputs:
% radi = radius
%
% Required inputs:
% bm   = beam structure
% 
% Optional inputs:
% 'norm'              : If set, indicates that the returned array
%                       contains the distances divided by the beam radius.
%                       This assumes the radius of the pilot tracer beam
%                       accurately reflects the size of the actual beam
%                       in the wavefront array, which will not be true
%                       in the case of significant aberrations.

% Revision history:
% 2005 Feb     jek  created idl routine
% 2016 Jun 21  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  [ny, nx] = size(bm.wf);
  ix2  = [-floor(nx / 2 ) : floor((nx - 1) / 2)].^2;
  iy2  = [-floor(ny / 2 ) : floor((ny - 1) / 2)].^2;
  [ax2, ay2]   = meshgrid(ix2, iy2);
  radi = sqrt(ax2 + ay2) * bm.dx;

  if (nargin > 1) & strcmp(norm, 'norm')
    radi = radi / prop_get_beamradius(bm);
  end

end                     % function prop_radius
