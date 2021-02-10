%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function mask = prop_rounded_rectangle(bm, rad, sx, sy, varargin)
%        mask = prop_rounded_rectangle(bm, rad, sx, sy, varargin)
% Return a 2-D array containing a rectangular mask (1 inside, 0 outside)
% with rounded corners.  This routine was created to allow modeling of
% the TPF-C promary mirror.
%
% Outputs:
% mask = 2-D array containing the mask
%
% Required inputs:
% bm   = wavefront structure (used for getting wavefront sampling info)
% rad  = radius of the rounded corners (which are 90 degree sections of
%        a circle) (m)
% sx   = size of the mask in x (m)
% sy   = size of the mask in y (m)
%
% Optional inputs:
% 'xc'                = rectangle center offset from the array center x (m)
% 'yc'                = rectangle center offset from the array center y (m)
%                       The default is a mask centered in the array

% Revision history:
% 2006 Feb     jek  created idl routine
% 2016 Jul 25  gmg  Matlab translation
% 2017 Apr 10  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  cx   = 0d0;                   % rectangle center X (m)
  cy   = 0d0;                   % rectangle center Y (m)

% Set values of internal parameters
  dx   = prop_get_sampling(bm); % spacing between array points x (m)
  dy   = prop_get_sampling(bm); % spacing between array points y (m)
  [ny, nx] = size(bm.wf);       % number of array points

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};  % rectangle center X (m)
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};  % rectangle center Y (m)
      otherwise
        error('prop_rounded_rectangle: Unknown keyword: %s\n', varargin{icav});
    end
  end

  radp = rad / dx;              % corner radius (pixels)
  cox  = sx / 2d0 - rad;        % circle offset x (m)
  coy  = sy / 2d0 - rad;        % circle offset y (m)

  mask = prop_rectangle(bm, sx, sy, 'cx', cx, 'cy', cy);

  circ = prop_ellipse(bm, rad, rad, 'cx', cx + cox, 'cy', cy + coy);
  oxp  = (cx + cox) / dx + fix(nx / 2) + 1;     % offset x (pixels)
  oyp  = (cy + coy) / dy + fix(ny / 2) + 1;     % offset y (pixels)
  x1   = fix(oxp);
  if x1 < 1
    x1   = 1;
  end
  x2   = fix(oxp + radp + 2);
  if x2 > nx
    x2   = nx;
  end
  y1   = fix(oyp);
  if y1 < 1
    y1   = 1;
  end
  y2   = fix(oyp + radp + 2);
  if y2 > nx
    y2   = nx;
  end
  mask(y1 : y2, x1 : x2) = 0d0;
  mask(y1 : y2, x1 : x2) = circ(y1 : y2, x1 : x2);

  circ = prop_ellipse(bm, rad, rad, 'cx', cx - cox, 'cy', cy + coy);
  oxp  = (cx - cox) / dx + fix(nx / 2) + 1;     % offset x (pixels)
  oyp  = (cy + coy) / dy + fix(ny / 2) + 1;     % offset y (pixels)
  x1   = fix(oxp - radp - 2);
  if x1 < 1
    x1   = 1;
  end
  x2   = fix(oxp);
  if x2 > nx
    x2   = nx;
  end
  y1   = fix(oyp);
  if y1 < 1
    y1   = 1;
  end
  y2   = fix(oyp + radp + 2);
  if y2 > nx
    y2   = nx;
  end
  mask(y1 : y2, x1 : x2) = 0d0;
  mask(y1 : y2, x1 : x2) = circ(y1 : y2, x1 : x2);

  circ = prop_ellipse(bm, rad, rad, 'cx', cx - cox, 'cy', cy - coy);
  oxp  = (cx - cox) / dx + fix(nx / 2) + 1;     % offset x (pixels)
  oyp  = (cy - coy) / dy + fix(ny / 2) + 1;     % offset y (pixels)
  x1   = fix(oxp - radp - 2);
  if x1 < 1
    x1   = 1;
  end
  x2   = fix(oxp);
  if x2 > nx
    x2   = nx;
  end
  y1   = fix(oyp - radp - 2);
  if y1 < 1
    y1   = 1;
  end
  y2   = fix(oyp);
  if y2 > nx
    y2   = nx;
  end
  mask(y1 : y2, x1 : x2) = 0d0;
  mask(y1 : y2, x1 : x2) = circ(y1 : y2, x1 : x2);

  circ = prop_ellipse(bm, rad, rad, 'cx', cx + cox, 'cy', cy - coy);
  oxp  = (cx + cox) / dx + fix(nx / 2) + 1;     % offset x (pixels)
  oyp  = (cy - coy) / dy + fix(ny / 2) + 1;     % offset y (pixels)
  x1   = fix(oxp);
  if x1 < 1
    x1   = 1;
  end
  x2   = fix(oxp + radp + 2);
  if x2 > nx
    x2   = nx;
  end
  y1   = fix(oyp - radp - 2);
  if y1 < 1
    y1   = 1;
  end
  y2   = fix(oyp);
  if y2 > nx
    y2   = nx;
  end
  mask(y1 : y2, x1 : x2) = 0d0;
  mask(y1 : y2, x1 : x2) = circ(y1 : y2, x1 : x2);

  end                           % prop_rounded_rectangle
