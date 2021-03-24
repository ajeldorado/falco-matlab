%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [bm, map] = prop_dm(bm, dmz0, dmcx, dmcy, spcg, varargin)
%        [bm, map] = prop_dm(bm, dmz0, dmcx, dmcy, spcg, varargin)
% Simulate a deformable mirror of specified actuator spacing,
% including the effects of the DM influence function.
%
% Outputs:
% bm   = beam structure
% map  = returns DM surface (not wavefront) on same grid as bm (m)
%
% Required inputs:
% bm   = beam structure
% dmz0 = either:
%        2D array containing the surface piston of each DM actuator (m)
%        or the name of a 2D FITS (Flexible Image Transport System)
%        image file containing the above
% dmcx, dmcy = the location of the optical axis (center of the
%        wavefront) on the DM in actuator units (0 to nAct - 1).
%        The center of the first actuator is (0.0, 0.0).
% spcg = defines the spacing between actuators (m);
%        must not be used when n_act_acroos_pupil is specified.
%
% Optional inputs:
% 'Fit'               : switch that tells routine that the values in dmz0
%                       are the desired surface heights rather than the
%                       commanded actuator heights, and so the routine
%                       should fit this map, accounting for actuator
%                       influence functions, to determine the necessary
%                       actuator heights.  An iterative error-minimizing
%                       loop is used for the fit.
% 'no_apply'          : if set, the DM pattern is not added to the
%                       wavefront.  Useful if the DM surface map is
%                       needed but should not be applied to the wavefront.
% 'n_act_across_pupil'= specifies the number of actuators that span the
%                       x-axis beam diameter.  If it is a whole number,
%                       the left edge of the left pixel is aligned with
%                       the left edge of the beam, and the right edge of
%                       the right pixel with the right edge of the beam.
%                       This determines the spacing and size of the
%                       actuators.  Should not be used when spcg value
%                       is specified.
% 'xtilt'             = specify the rotation of the DM surface with
% 'ytilt'               respect to the wavefront plane, with the origin
% 'ztilt'               at the center of the wavefront.  The DM surface
%                       is interpolated and orthogonally projected onto
%                       the wavefront grid.  The coordinate system
%                       assumes that the wavefront and initial DM
%                       surface are in the X,Y plane with a lower left
%                       origin with Z towards the observer.  The
%                       rotations are left handed.  The default rotation
%                       order is X, Y, then Z unless the zyx switch is
%                       set. (default: 0, 0, 0)
% 'zyx'               : specifies the rotation order if two or more of
%                       xtilt, ytilt, or ztilt are specified.
%                       The default is X, Y, then Z rotations.
%                       Set for for Z, Y, then X rotations.

% 2005 Feb     jek  created idl routine 
% 2007 Jun     jek  fixed bug in interpolation of smoothed DM surface &
%                   changed default from smoothed to non-smoothed
% 2014 Mar     jek  implemented faster method of creating DM surface
% 2014 Jun 26  gmg  Matlab translation
% 2015 Oct     jek  Implemented DM surface tilting, rotation
% 2015 Dec 04  gmg  Matlab translation
% 2017 Feb 15  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ischar(dmz0)               % then open 2D FITS image file
    dmz  = fitsread(dmz0);
  else
    dmz  = dmz0;
  end

  Fit  = 0;             % values in dmz0 are commanded actuator heights
  nAct = 0;             % number of actuators across pupil
  NoAp = 0;             % the DM pattern is added to the wavefront
  tlt  = zeros(1, 3);   % set default tilts to 0 degrees
  zyx  = 0;             % rotation order is X, Y, then Z rotations

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'fit'}
        Fit  = 1;       % values in dmz0 are desired surface heights
      case {'noap', 'no_apply'}
        NoAp = 1;       % the DM pattern is not added to the wavefront
      case {'nact', 'n_act_across_pupil'}
        icav = icav + 1;
        nAct = varargin{icav};  % number of actuators across pupil
      case {'tltx', 'xtilt'}
        icav = icav + 1;
        tlt(1) = varargin{icav};% X rotation of the DM surface (deg)
      case {'tlty', 'ytilt'}
        icav = icav + 1;
        tlt(2) = varargin{icav};% Y rotation of the DM surface (deg)
      case {'tltz', 'ztilt'}
        icav = icav + 1;
        tlt(3) = varargin{icav};% Z rotation of the DM surface (deg)
      case {'xyz'}
        zyx  = 0;       % rotation order is X, Y, then Z rotations
      case {'zyx'}
        zyx  = 1;       % rotation order is Z, Y, then X rotations
      otherwise
        error('prop_dm: Unknown keyword: %s\n', varargin{icav});
    end
  end

% Read the influence function data from a FITS file
  info = fitsinfo('influence_dm5v2.fits');
  inf  = fitsread('influence_dm5v2.fits');

  [ldef, idef] = ismember('NAXIS1' , info.PrimaryData.Keywords(:, 1));
  ifnx = info.PrimaryData.Keywords{idef, 2};    % inf func number of pixels x

  [ldef, idef] = ismember('NAXIS2' , info.PrimaryData.Keywords(:, 1));
  ifny = info.PrimaryData.Keywords{idef, 2};    % inf func number of pixels y

  [ldef, idef] = ismember('P2PDX_M', info.PrimaryData.Keywords(:, 1));
  ifdx = info.PrimaryData.Keywords{idef, 2};    % inf func spacing x (m)

  [ldef, idef] = ismember('P2PDY_M', info.PrimaryData.Keywords(:, 1));
  ifdy = info.PrimaryData.Keywords{idef, 2};    % inf func spacing y (m)

  [ldef, idef] = ismember('C2CDX_M', info.PrimaryData.Keywords(:, 1));
  acdx = info.PrimaryData.Keywords{idef, 2};    % actuator spacing x (m)

  [ldef, idef] = ismember('C2CDY_M', info.PrimaryData.Keywords(:, 1));
  acdy = info.PrimaryData.Keywords{idef, 2};    % actuator spacing y (m)

  ifcx = floor(ifnx / 2) + 1;           % inf func center pixel x
  ifcy = floor(ifny / 2) + 1;           % inf func center pixel y
  raix = round(acdx / ifdx);            % ratio actuator / inf func dx
  raiy = round(acdy / ifdy);            % ratio actuator / inf func dy

  if spcg ~=  0.0 & nAct ~=  0
    error('Proper:PROP_DM', ...
          'User cannot specify both actuator spacing and number of actuators across pupil.\n');
  end

  if (spcg ==  0.0 & nAct ==  0)
    error('Proper:PROP_DM', ...
          'User must specify either actuator spacing or number of actuators across pupil.\n');
  end

% Set the real DM actuator spacing
  if nAct ~= 0
    actd = 2.0 * prop_get_beamradius(bm) / nAct;
  else
    actd = spcg;
  end

% Scale the influence function sampling to the specified DM actuator spacing
  ifdx = ifdx * actd / acdx;            % actual inf func spacing x (m)
  ifdy = ifdy * actd / acdy;            % actual inf func spacing y (m)

  if Fit == 1
% Then calculate the actuator positions to fit the given DM surface
% Create the meshgrids for the influence function inf
    [cfx, cfy] = meshgrid([ 1 : ifnx], [ 1 : ifny]);
% Create the meshgrids for the interpolated influence function infk
    [ckx, cky] = meshgrid([-2 :  2  ], [-2 :  2  ]);
    ckx  = ckx * actd;
    cky  = cky * actd;
% Calculate the interpolated influence function infk
    infk = interp2(cfx, cfy, inf, ckx/ifdx + ifcx, cky/ifdy + ifcy, 'cubic');
% Calculate the actuator positions to fit the given DM surface
    dmcz = prop_fit_dm(dmz, infk);
  else
% Use the given actuator positions
    dmcz = dmz;
  end

  [acny, acnx] = size(dmz);             % number of actuators in y and x
% Create subsampled DM grid
  mrgx = raix * 9;                      % margin on each end of each row
  mrgy = raiy * 9;                      % margin on each end of each column
  dmnx = raix * acnx + mrgx * 2;        % number of DM grid points x
  dmny = raiy * acny + mrgy * 2;        % number of DM grid points y
  lgx1 = fix(raix/2) + mrgx + 1;        % index of  1st actuator center x
  lgy1 = fix(raiy/2) + mrgy + 1;        % index of  1st actuator center y
  lgx2 = lgx1 + raix * (acnx-1);        % index of last actuator center x
  lgy2 = lgy1 + raiy * (acny-1);        % index of last actuator center y
  dmg  = zeros(dmny, dmnx);             % initialize subsampled DM grid
% Fill the subsampled DM grid with the actuator amplitudes
  dmg(lgy1 : raiy : lgy2, lgx1 : raix : lgx2) = dmcz;
% Do the convolution of the actuator amplitudes with the influence function
  dmg  = conv2(dmg, inf, 'same');

  [ny, nx] = size(bm.wf);       % number of points in wavefront array

% 3D rotate DM grid and project orthogonally onto wavefront

% Calculate grid dimensions (pix) projected onto wavefront
  xdim = round(sqrt(2.0) * dmnx * ifdx / bm.dx);
  if xdim > nx
    xdim = nx;
  end
  xd2  = fix(xdim / 2) + 1;
  ydim = round(sqrt(2.0) * dmny * ifdy / bm.dx);
  if ydim > ny
    ydim = ny;
  end
  yd2  = fix(ydim / 2) + 1;

  cx   = ([1 : xdim] - xd2) * bm.dx;
  cy   = ([1 : ydim] - yd2) * bm.dx;
  [cxm, cym] = meshgrid(cx, cy);

  sa   = sind(tlt(1));
  ca   = cosd(tlt(1));
  sb   = sind(tlt(2));
  cb   = cosd(tlt(2));
  sg   = sind(tlt(3));
  cg   = cosd(tlt(3));

  if zyx == 0
    rotm = [               cb * cg,               -cb * sg,       sb, 0.0; ...
            ca * sg + sa * sb * cg, ca * cg - sa * sb * sg, -sa * cb, 0.0; ...
            sa * sg - ca * sb * cg, sa * cg + ca * sb * sg,  ca * cb, 0.0; ...
                               0.0,                    0.0,      0.0, 1.0];
  else
    rotm = [ cb * cg, sa * sb * cg - ca * sg, ca * sb * cg + sa * sg, 0.0; ...
             cb * sg, sa * sb * sg + ca * cg, ca * sb * sg - sa * cg, 0.0; ...
            -sb,      sa * cb,                ca * cb,                0.0; ...
                 0.0,                    0.0,                    0.0, 1.0];
  end

% Forward project a square
  edge = [-1.0, -1.0,  0.0,  0.0; ...
           1.0, -1.0,  0.0,  0.0; ...
           1.0,  1.0,  0.0,  0.0; ...
          -1.0,  1.0,  0.0,  0.0];
  xyzn = edge * rotm;   % had to reverse matrix mult. order to match IDL

% Determine backwards projection for screen-raster-to-DM-surface computation
% Had to reverse and increment indices to match IDL
  dxdx = (xyzn(1, 1) - xyzn(2, 1)) / (edge(1, 1) - edge(2, 1));
  dxdy = (xyzn(2, 1) - xyzn(3, 1)) / (edge(2, 2) - edge(3, 2));
  dydx = (xyzn(1, 2) - xyzn(2, 2)) / (edge(1, 1) - edge(2, 1));
  dydy = (xyzn(2, 2) - xyzn(3, 2)) / (edge(2, 2) - edge(3, 2));

  xs   = (cxm/dxdx - cym*dxdy / (dxdx*dydy)) / (1.0 - dydx*dxdy / (dxdx*dydy));
  ys   = (cym/dydy - cxm*dydx / (dxdx*dydy)) / (1.0 - dydx*dxdy / (dxdx*dydy));

  xdm  = (xs + dmcx * actd) / ifdx + lgx1;
  ydm  = (ys + dmcy * actd) / ifdy + lgy1;

% Calculate the interpolated DM grid (set extrapolated values to 0.0)
  grid = interp2(dmg, xdm, ydm, 'cubic', 0.0);

% Create the limits for the map
  icx  = fix(nx / 2) + 1;       % wf array center index x
  icy  = fix(ny / 2) + 1;       % wf array center index y
  isx1 = icx  - xd2  + 1;       %  1st point index in wf subgrid x
  isy1 = icy  - yd2  + 1;       %  1st point index in wf subgrid y
  isx2 = isx1 + xdim - 1;       % last point index in wf subgrid x
  isy2 = isy1 + ydim - 1;       % last point index in wf subgrid y
  igx1 = 1;                     %  1st point index in DM grid x
  igy1 = 1;                     %  1st point index in DM grid y
  igx2 = xdim;                  % last point index in DM grid x
  igy2 = ydim;                  % last point index in DM grid y
  if isx1 < 1           % clip interpolated DM grid at wf array bounds
    igx1 = 1 + (1 - isx1);
    isx1 = 1;
  end
  if isy1 < 1           % clip interpolated DM grid at wf array bounds
    igy1 = 1 + (1 - isy1);
    isy1 = 1;
  end
  if isx2 > nx          % clip interpolated DM grid at wf array bounds
    igx2 = xdim - (isx2 - nx);
    isx2 = nx;
  end
  if isy2 > ny          % clip interpolated DM grid at wf array bounds
    igy2 = ydim - (isy2 - ny);
    isy2 = ny;
  end

  map  = zeros(ny, nx);
  map(isy1 : isy2, isx1 : isx2) = grid(igy1 : igy2, igx1 : igx2);

  if NoAp == 0
    bm   = prop_add_phase(bm, 2.0 * map);
  end
end                     % function prop_dm
