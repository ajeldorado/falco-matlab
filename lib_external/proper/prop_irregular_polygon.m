%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function apm = prop_irregular_polygon(bm, vx, vy, varargin)
%        apm = prop_irregular_polygon(bm, vx, vy, varargin)
% Return an image containing a filled (interior value is 1.0) convex
% polygon with antialiased edges.  If the dark switch is set, then the
% interior is set to 0.0.  NOTE: This routine only works on convex
% polygons and will produce incorrect results given anything else.
%
% Outputs:
% apm  = aperture mask containing antialiased filled polygon
%
% Required inputs:
% bm   = beam structure (used only to obtain the beam radius and sampling
%        of the current wavefront; the wavefront array is not modified)
% vx,vy= vector arrays of the same dimension containing the coordinates
%        of the convex polygon vertices (meters unless norm = 1, then
%        fraction of beam diameter)  The center of the wavefront array
%        is designated to have the coordinates (0, 0).  The vertices
%        must be in sequential order, either clockwise or counter-
%        clockwise.  The coordinates of the first and last vertices
%        should be the same (i.e., close the polygon).
%
% Optional inputs:
% 'dark'              : draw a dark poly. (0 inside, 1 outside)
%                       (default is opposite way)
% 'norm'              : vertex coordinates are normalized to beam radius

% 2005 Feb     jek  created idl routine
% 2016 Aug 12  gmg  Matlab translation
% 2017 Mar 21  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

% Set default values of input parameters
  dark = 0;             % no dark center
  norm = 0;             % vertex coordinates in meters

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'dark'}
        dark = 1;               % dark center
      case {'norm'}
        norm = 1;               % vertex coordinates normalized to beam radius
      otherwise
        error('prop_irregular_polygon: Unknown keyword: %s\n', varargin{icav});
    end
  end

  dx   = prop_get_sampling(bm); % spacing between points in x (m)
  dy   = prop_get_sampling(bm); % spacing between points in y (m)
  magn = double(antialias_subsampling);                    % subsampling factor at edges; must be odd
  [ny, nx] = size(bm.wf);               % number of pixels
  prx  = prop_get_beamradius(bm) / dx;  % beam radius x in pixels
  pry  = prop_get_beamradius(bm) / dy;  % beam radius y in pixels

  if size(vx, 1) > 1
    nvrt = size(vx, 1);         % number of vertices
  else
    nvrt = size(vx, 2);         % number of vertices
  end

% Drawing algorithm needs vertices to be clockwise-ordered.
% Test for clockwise-ness of vertices; for a convex polygon, if the
% cross product of adjacent edges is positive, it is counter-clockwise.
  cp   = (vx(2)-vx(1)) * (vy(3)-vy(2)) - (vy(2)-vy(1)) * (vx(3)-vx(2));
  if cp > 0
    vx   = vx([nvrt : -1 : 1]);
    vy   = vy([nvrt : -1 : 1]);
  end

  if norm == 1
    vx   = vx * prx + fix(nx / 2) + 1;
    vy   = vy * pry + fix(ny / 2) + 1;
  else
    vx   = vx /  dx + fix(nx / 2) + 1;
    vy   = vy /  dy + fix(ny / 2) + 1;
  end

  apm  = zeros(ny, nx);                 % initialize mask to zeros

  left = find(vy == min(vy));
  left = left(find(vx(left) == min(vx(left))));
  left = left(1);                       % index of first left corner
  if left ~= nvrt
    ltnt = left + 1;                    % index of next left corner
  else
    ltnt = 1;                           % index of next left corner
  end
  rght = left;                          % index of first right corner
  if rght ~= 1
    rtnt = rght - 1;                    % index of next right corner
  else
    rtnt = nvrt;                        % index of next right corner
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
        if xltp >= 1
          apm(ipy, xltp) = apm(ipy, xltp) + magn * ((xltp + 0.5) - xlt);
        end
        if xrtp <= nx
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
end                     % function prop_irregular_polygon
