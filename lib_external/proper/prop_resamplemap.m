%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [mapo] = prop_resamplemap(bm, mapi, pxsc, cpx, cpy, xshift, yshift)
%        [mapo] = prop_resamplemap(bm, mapi, pxsc, cpx, cpy, xshift, yshift)
% Interpolate input map using cubic convolution onto grid with same size
% and sampling as the current wavefront array.  Optionally shift the map.
%
% Outputs:
% mapo = output map
%
% Required inputs:
% bm   = beam structure
% mapi = aberration map to be resampled
% pxsc = spacing of "mapi" (m)
% cpx  = pixel coordinate of mapi center X (0, 0 is center of 1st pixel)
% cpy  = pixel coordinate of mapi center Y (0, 0 is center of 1st pixel)
%
% Optional inputs:
% xshift              = amount to shift map X (m)
% yshift              = amount to shift map Y (m)
%
% Intended for internal use only by prop_* routines.

% 2005 Feb     jek  created idl routine
% 2014 Jul 29  gmg  Matlab translation
% 2017 Feb 23  gmg  Revised for keyword/value for optional inputs
% 2017 Jun 19  gmg  Revised to replace keyword/value for optional inputs
%                   with ordered optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if nargin < 6
    xshift = 0.0;               % set default map shift X to 0.0
  end
  if nargin < 7
    yshift = 0.0;               % set default map shift Y to 0.0
  end

  [ny, nx] = size(bm.wf);       % number of points in wavefront array
  o1x  = ([1 : nx] - fix(nx / 2) - 1)  * bm.dx / pxsc;
  o1x  = o1x + cpx - xshift / pxsc;
  o1y  = ([1 : ny] - fix(ny / 2) - 1)' * bm.dx / pxsc;
  o1y  = o1y + cpy - yshift / pxsc;

  mapo = interp2(mapi, o1x, o1y, 'cubic', 0.0);
end                     % function prop_resamplemap
