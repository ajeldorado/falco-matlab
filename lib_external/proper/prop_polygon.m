%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function img = prop_polygon(bm, nvrt, rad, varargin)
%        img = prop_polygon(bm, nvrt, rad, varargin)
% Return an image containing a filled, symmetric (vertices are the same
% distance from the center) polygon with anti-aliased edges.  A clear
% polygon is drawn (one inside, zero outside), unless "dark" is set.
% Unless "rot" is defined, the polygon is drawn with one vertex
% positioned along the +X axis. relative to the polygon center.
%
% Outputs:
% img  = image containing the polygon
%
% Required inputs:
% bm   = Current beam structure.  This is used only to obtain the
%        sampling of the current wavefront; the wavefront array is not
%        modified.
% nvrt = number of vertices; must be 3 or greater
% rad  = distance from any vertex to the center of the polygon (same
%        for all vertices).  By default this is in meters.  If "norm" is
%        set, then this is relative to the current beam radius.
%
% Optional inputs:
% 'xc'                = polygon center X relative to the wavefront center X.
% 'yc'                = polygon center Y relative to the wavefront center Y.
%                       These are in meters unless "norm" is set, in
%                       which case they are relative to the current beam
%                       radius.  The default is for the polygon to be
%                       centered in the array.
% 'rotation'          = counter-clockwise angle to rotate the polygon
%                       (degrees).  By default, one vertex lies along
%                       the +X axis relative to the polygon center.
% 'dark'              : if set, then the interior of the polygon will be
%                       set to 0.0 and the exterior to 1.0.
% 'norm'              : if set, indicates that the radius and polygon
%                       center values are relative to the current beam radius.

% 2005 Feb     jek  created idl routine
% 2016 Jul 27  gmg  Matlab translation
% 2017 Mar 16  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

% Set default values of input parameters
  cx   = 0.0;           % pixel coordinates of polygon center X
  cy   = 0.0;           % pixel coordinates of polygon center Y
  dark = 0;             % interior of polygon = 1, exterior = 0
  norm = 0;             % radius and center are in m
  rot  = 0.0;           % counter-clockwise rotation of polygon (deg)

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};  % pixel coordinates of polygon center X
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};  % pixel coordinates of polygon center Y
      case {'dark'}
        dark = 1;               % dark center
      case {'norm'}
        norm = 1;               % radius and center relative to beam radius
      case {'rot', 'rotation'}
        icav = icav + 1;
        rot  = varargin{icav};  % counter-clockwise rotation of polygon (deg)
      otherwise
        error('prop_polygon: Unknown keyword: %s\n', varargin{icav});
    end
  end

  dx   = prop_get_sampling(bm);         % spacing between pixels x (m)
  dy   = prop_get_sampling(bm);         % spacing between pixels y (m)
  brxp = prop_get_beamradius(bm) / dx;  % beam radius x in pixels
  bryp = prop_get_beamradius(bm) / dy;  % beam radius x in pixels
  [ny, nx] = size(bm.wf);               % size of wavefront array

  cxp  = fix(nx / 2) + 1;               % polygon center x (pixels)
  cyp  = fix(ny / 2) + 1;               % polygon center y (pixels)
  radx = rad / dx;                      % polygon radius x (pixels)
  rady = rad / dy;                      % polygon radius y (pixels)
  if norm == 1
    cxp  = cxp + cx * brxp;
    cyp  = cyp + cy * bryp;
    radx = rad * brxp;
    rady = rad * bryp;
  else
    cxp  = cxp + cx / dx;
    cyp  = cyp + cy / dy;
  end

  magn = double(antialias_subsampling);            % subsampling factor at edges; must be odd
  ang  = -360d0 * [0 : nvrt - 1] / nvrt;
  vx0  = cosd(ang);             % polygon vertex coordinate x
  vy0  = sind(ang);             % polygon vertex coordinate y
% Calculate vertex coordinates of rotated and displaced polygon
  vx   = cxp + radx * (vx0 * cosd(rot) - vy0 * sind(rot));
  vy   = cyp + rady * (vx0 * sind(rot) + vy0 * cosd(rot));
  img  = zeros(ny, nx);

  left = find(vy == min(vy));
  left = left(find(vx(left) == min(vx(left))));
  left = left(1);               % index of lower left corner (1 - nvrt)
  if left ~= nvrt
    ltnt = left + 1;            % index of next left corner (1 - nvrt)
  else
    ltnt = 1;                   % index of next left corner (1 - nvrt)
  end
  rght = left;
  if rght ~= 1
    rtnt = rght - 1;            % index of next right corner (1 - nvrt)
  else
    rtnt = nvrt;                % index of next right corner (1 - nvrt)
  end

  ivyn = round(min(vy));        % index of vertex minimum y
  if ivyn < 1
    ivyn = 1;
  end
  ivyx = round(max(vy));        % index of vertex maximum y
  if ivyx > ny
    ivyx = ny;
  end

  for ipy = ivyn : ivyx         % index of pixel y
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
          img(ipy, xltp) = img(ipy, xltp) + magn * ((xltp + 0.5) - xlt);
        end
        if xrtp <= nx
          img(ipy, xrtp) = img(ipy, xrtp) + magn * (xrt - (xrtp - 0.5));
        end
        if (xrtp - xltp) > 1
          xlp1 = max(xltp + 1,  1);
          xrp1 = min(xrtp - 1, nx);
          img(ipy, xlp1 : xrp1) = img(ipy, xlp1 : xrp1) + magn;
        end
      else
        if ((xltp >= 1) & (xltp <= nx))
          img(ipy, xltp) = img(ipy, xltp) + magn * (xrt - xlt);
        end
      end                       % if xltp ~= xrtp
    end                         % for suby = 0 : (magn - 1)
  end                           % for ipy = ivyn : ivyx

  img  = img / magn^2;
  if dark == 1
    img  = 1.0 - img;
  end
end                     % prop_polygon
