%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function apm = prop_rectangle(bm, sx, sy, varargin)
%        apm = prop_rectangle(bm, sx, sy, varargin)
% Return an image containing an antialiased, filled rectangle
%
% Outputs:
% apm  = aperture mask containing antialiased filled rectangle
%
% Required inputs:
% bm   = beam structure
% sx   = size along x (m unless 'norm', then fraction of beam radius)
% sy   = size along y (m unless 'norm', then fraction of beam radius)
%
% Optional inputs:
% 'xc'                = center of rectangle relative to wf center X
%                       (m unless 'norm')
% 'yc'                = center of rectangle relative to wf center Y
%                       (m unless 'norm')
% 'rotation'          = counter-clockwise rotation of rectangle about
%                       its center (deg)
% 'dark'              : draw a dark rectangle (0 inside, 1 outside)
%                       (default is opposite way)
% 'norm'              = sizes and center coordinates are normalized
%                       to beam radius

% 2005 Feb     jek  created idl routine
% 2015 Nov 11  gmg  Matlab translation
% 2017 Feb 23  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

% Set default values of input parameters
  cx   = 0.0;                           % center of rectangle X
  cy   = 0.0;                           % center of rectangle Y
  dark = 0;
  norm = 0;
  rot  = 0.0;                           % rotation (degrees)

% Set values of internal parameters
  dx   = prop_get_sampling(bm);         % spacing between points in x (m)
  dy   = prop_get_sampling(bm);         % spacing between points in y (m)
  magn = double(antialias_subsampling);         % subsampling factor at edges; must be odd
  [ny, nx] = size(bm.wf);       % number of pixels in wavefront array
  prx  = prop_get_beamradius(bm) / dx;  % beam radius x in pixels
  pry  = prop_get_beamradius(bm) / dy;  % beam radius y in pixels

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};
      case {'dark'}
        dark = 1;
      case {'norm'}
        norm = 1;
      case {'rot', 'rotation'}
        icav = icav + 1;
        rot  = varargin{icav};  % rotation (degrees)
      otherwise
        error('prop_rectangle: Unknown keyword: %s\n', varargin{icav});
    end
  end

  cpx  = floor(nx / 2 + 1) + cx / dx;   % center X (pixels)
  cpy  = floor(ny / 2 + 1) + cy / dy;   % center Y (pixels)
  radx = 0.5 * sx / dx;                 % radius in X (pixels)
  rady = 0.5 * sy / dy;                 % radius in Y (pixels)

  if norm == 1
    cpx  = floor(nx / 2 + 1) + cx * prx;
    cpy  = floor(ny / 2 + 1) + cy * pry;
    radx = 0.5 * sx * prx;              % radius in X (pixels)
    rady = 0.5 * sy * pry;              % radius in Y (pixels)
  end

% Calculate rectangle vertex x and y coordinates before rotation and shift
  vx0  = [-radx, -radx,  radx,  radx];  % (pixels)
  vy0  = [-rady,  rady,  rady, -rady];  % (pixels)
  nvrt = 4;                             % number of vertices

% Calculate rectangle vertex x and y coordinates after rotation and shift
  vx   = vx0 * cosd(rot) - vy0 * sind(rot) + cpx;       % (pixels)
  vy   = vx0 * sind(rot) + vy0 * cosd(rot) + cpy;       % (pixels)

  apm  = zeros(ny, nx);                 % initialize mask to zeros

  left = find(vy == min(vy));
  left = left(find(vx(left) == min(vx(left))));
  left = left(1);                       % index of lower left corner (1 - 4)
  if left ~= nvrt
    ltnt = left + 1;                    % index of next left corner (1 - 4)
  else
    ltnt = 1;                           % index of next left corner (1 - 4)
  end
  rght = left;
  if rght ~= 1
    rtnt = rght - 1;                    % index of next right corner (1 - 4)
  else
    rtnt = nvrt;                        % index of next right corner (1 - 4)
  end

  ivyn = round(min(vy));                % index of vertex minimum y
  if ivyn < 1
    ivyn = 1;
  end
  ivyx = round(max(vy));                % index of vertex maximum y
  if ivyx > ny
    ivyx = ny;
  end

  for ipy = ivyn : ivyx                 % index of pixel y
    for suby = 0 : (magn - 1)
% Calculate y coordinate of current edge pixel (pixels)
      py   = ipy - 0.5 + (0.5 + suby) / magn;

      if py < vy(left)
        continue;
      end
      if py > max(vy)
        break;
      end
      if py >= vy(ltnt)
        left = ltnt;
        if left ~= nvrt
          ltnt = left + 1;
        else
          ltnt = 1;
        end
      end

      if py >= vy(rtnt)
        rght = rtnt;
        if rght ~= 1
          rtnt = rght - 1;
        else
          rtnt = nvrt;
        end
      end

      dylt = vy(ltnt) - vy(left);
      if dylt ~= 0
        dxlt = vx(ltnt) - vx(left);
        xlt  = dxlt / dylt * (py - vy(left)) + vx(left);
      else
        xlt  = vx(left);
      end

      dyrt = vy(rtnt) - vy(rght);
      if dyrt ~= 0
        dxrt = vx(rtnt) - vx(rght);
        xrt  = dxrt / dyrt * (py - vy(rght)) + vx(rght);
      else
        xrt  = vx(rght);
      end

      xltp = round(xlt);
      xrtp = round(xrt);

      if xltp ~= xrtp
        if ((xltp >= 1) & (xltp <= nx))
          apm(ipy, xltp) = apm(ipy, xltp) + magn * ((xltp + 0.5) - xlt);
        end
        if ((xrtp >= 1) & (xrtp <= nx))
          apm(ipy, xrtp) = apm(ipy, xrtp) + magn * (xrt - (xrtp - 0.5));
        end
        if (xrtp - xltp) > 1
          xlp1 = max(xltp + 1,  1);
          xrp1 = min(xrtp - 1, nx);
          apm(ipy, xlp1 : xrp1) = apm(ipy, xlp1 : xrp1) + magn;
        end
      else
        if ((xltp >= 1) & (xltp <= nx))
          apm(ipy, xltp) = apm(ipy, xltp) + magn * (xrt - xlt);
        end
      end                               % if xltp ~= xrtp
    end                                 % for suby = 0 : (magn - 1)
  end                                   % for ipy = ivyn : ivyx

  apm  = apm / magn^2;
  if dark == 1
    apm  = 1.0 - apm;
  end
