%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function arro = prop_rotate(arri, angd, varargin)
%        arro = prop_rotate(arri, angd, varargin)
% Rotate and shift an array via interpolation (bilinear by default).
% Returns a rotated and shifted array with the same dimensions as the
% input array.
%
% Outputs:
% arro = 2D output array
%
% Required inputs:
% arri = array to be rotated
% angd = angle to rotate array counter-clockwise (degrees)
%
% Optional inputs:
% 'meth'              = interpolation method
%         = 'nearest' = nearest neighbor interpolation
%         = 'linear'  = bilinear interpolation (default)
%         = 'spline'  = spline interpolation
%         = 'cubic'   = bicubic interpolation as long as the data is
%                       uniformly spaced, otherwise the same as 'spline'
%                     = same as IDL interpolate routine with CUBIC = -.5
% 'missing'           = value to set extrapolated pixels (default = 0.0)
% 'xshift'            = amount to shift arro in x direction (pixels)
%                       (default = 0.0)
% 'yshift'            = amount to shift arro in y direction (pixels)
%                       (default = 0.0)
% 'xc'                = pixel coordinates of arri center
%                       (default = fix(nx / 2) + 1)
% 'yc'                = pixel coordinates of arri center
%                       (default = fix(ny / 2) + 1)

% 2005 Feb     jek  created idl routine
% 2014 Aug 18  gmg  Matlab translation; works correctly for non-square arrays
% 2017 Mar 10  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  [ny, nx] = size(arri);
  cx = fix(nx / 2) + 1;
  cy = fix(ny / 2) + 1;
  extr = 0.0;
  meth = 'linear';
  sx   = 0.0;
  sy   = 0.0;

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};  % input array center X (pixels)
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};  % input array center Y (pixels)
      case {'extr', 'missing'}
        icav = icav + 1;
        extr = varargin{icav};  % extrapolated pixel value
      case {'cubic'}
        meth = 'cubic';         % interpolation method
      case {'meth'}
        icav = icav + 1;
        meth = varargin{icav};  % interpolation method
      case {'sx', 'xshift'}
        icav = icav + 1;
        sx   = varargin{icav};  % output array shift X (pixels)
      case {'sy', 'yshift'}
        icav = icav + 1;
        sy   = varargin{icav};  % output array shift Y (pixels)
      otherwise
        error('prop_rotate: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  if nargin < 2
    error('Proper:PROP_ROTATE', ...
          'Must specify input array and rotation angle.\n');
  end

  [aox, aoy] = meshgrid([1 : nx], [1 : ny]);    % output array coordinates
% Minus sign for angle is needed for compatibility with IDL prop_rotate.pro
  cosa = cosd(-angd);
  sina = sind(-angd);
% Transform output array coordinates to input array coordinates
  aix  = (aox - cx - sx) * cosa - (aoy - cy - sy) * sina + cx;
  aiy  = (aox - cx - sx) * sina + (aoy - cy - sy) * cosa + cy;
  arro = interp2(arri, aix, aiy, meth, extr);
end                     % function prop_rotate
