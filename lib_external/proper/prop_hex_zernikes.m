%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function abm = prop_hex_zernikes(izp, zc, nx, dx, hr, varargin)
%        abm = prop_hex_zernikes(izp, zc, nx, dx, hr, varargin)
% Return a 2D array that is the sum of specified hexagonal Zernike
% polynomials.
% Note: The aberration map extends beyond the hexagonal region inside
% of which it is properly normalized.  If adding together a number of
% aberration maps generated with this routine, to create a wavefront map
% for a hexagonally segmented aperture for example, then each segment's
% map must be multiplied by its corresponding aperture prior to addition.
% prop_hex_wavefront will do all this.  The normalization error
% (difference between the requested and actual RMS % of the aberration)
% is about 0.5%, except for Z18 (1.5%) and Z22 (1%).
%
% Explanation: The user provides an array of Zernike polynomial indices
% (1 to 22) and an array of corresponding polynomial coefficients that
% specify the RMS amount of each aberration.  The polynomials are
% circular but normalized over a hexagonal aperture with the zero-angle
% axis passing through a vertex of the hexagon.  The order of the
% Zernikes is the same as the first 22 circular Zernikes used in PROPER
% (see prop_print_zernikes or prop_zernikes).
% Derivation of these polynomials is provided by
% V.N. Mahajan & G.-m. Dai, J. Opt. Soc. Am. A, 24, 9, 2994-3016 (2007)
%
% Outputs:
% abm  = the 2D array of summed Zernike polynomials (meters)
%
% Required inputs:
% izp  = array of Zernike polynomial indices (1 to 22)
% zc   = array of coefficients giving the RMS amount of aberration for
%        the corresponding Zernike polynomials specified by the indices
%        in "izp". (meters)
% nx   = size of array grid (nx by nx pixels)
% dx   = spacing of array grid (meters)
% hr   = distance from the center of the hexagon to a vertex (meters)
%
% Optional inputs:
% 'rotation'          = if given, the aberration array will be rotated
%                       by "ang" degrees about (hx, hy) in the clockwise
%                       direction (first pixel in lower left). (degrees)
% 'xhex'              = offset of the center of the hexagon from the
% 'yhex'                center of the wavefront.  If not specified, the
%                       hexagon is assumed to be centered in the
%                       wavefront at (hx, hy) = (0, 0). (meters)

% Revision history:
% 2005 Feb     jek
% 2016 Jul 11  gmg  Translated to Matlab
% 2017 Mar 16  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  ang  = 0.0;           % clockwise rotation of aberration array (deg)
  hx   = 0.0;           % hexagon center to wavefront center X (m)
  hy   = 0.0;           % hexagon center to wavefront center Y (m)

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'hx', 'xhex'}
        icav = icav + 1;
        hx   = varargin{icav};  % hexagon center to wavefront center X (m)
      case {'hy', 'yhex'}
        icav = icav + 1;
        hy   = varargin{icav};  % hexagon center to wavefront center Y (m)
      case {'ang', 'rotation'}
        icav = icav + 1;
        ang  = varargin{icav};  % clockwise rotation of aberration array (deg)
      otherwise
        error('prop_hex_zernikes: Unknown keyword: %s\n', varargin{icav});
    end
  end

  dy   = dx;                    % spacing of array grid in y (meters)
  ny   = nx;                    % number of points in output array y
  iacx = fix(nx / 2) + 1;       % x index of center of output array
  iacy = fix(ny / 2) + 1;       % y index of center of output array

  hxp  = round(hx / dx);        % center of hexagon x (pixels)
  hyp  = round(hy / dy);        % center of hexagon y (pixels)
  hrxp = fix(hr / dx) + 2;      % radius of hexagon x (pixels)
  hryp = fix(hr / dy) + 2;      % radius of hexagon y (pixels)

  izx1 = hxp  - hrxp;           % x index of minimum point in zer
  iax1 = iacx + izx1;           % x index of minimum point in abm
  if iax1 <  1
    izx1 = izx1 - iax1 +  1;
    iax1 =  1;
  end
    
  izx2 = hxp  + hrxp;           % x index of maximum point in zer
  iax2 = iacx + izx2;           % x index of maximum point in abm
  if iax2 > nx
    izx2 = izx2 - iax2 + nx;
    iax2 = nx;
  end

  izy1 = hyp  - hryp;           % y index of minimum point in zer
  iay1 = iacy + izy1;           % y index of minimum point in abm
  if iay1 <  1
    izy1 = izy1 - iay1 +  1;
    iay1 =  1;
  end
    
  izy2 = hyp  + hryp;           % y index of maximum point in zer
  iay2 = iacy + izy2;           % y index of maximum point in abm
  if iay2 > ny
    izy2 = izy2 - iay2 + ny;
    iay2 = ny;
  end

  izxs = izx2 - izx1 + 1;       % zer array size x (pixels)
  izys = izy2 - izy1 + 1;       % zer array size y (pixels)

  [cx, cy] = meshgrid((izx1 : izx2), (izy1 : izy2));
  cx   = cx * dx - hx;          % x coordinates of grid points (m)
  cy   = cy * dy - hy;          % y coordinates of grid points (m)

  cr   = sqrt(cx.^2 + cy.^2) / hr;          % normalized radius of grid points
  ct   = atan2(cy, cx) - ang * pi / 180d0;  % azimuth angle of grid points (rad)

  zer  = zeros(izys, izxs);     % 2D array of Zernike values in hexagon (m)

  cz   = zeros(22, 4);          % coefficients in Zernike polynomials

  cz( 2, 1) =    2d0 * sqrt( 6d0 /        5d0);
  cz( 3, 1) =    2d0 * sqrt( 6d0 /        5d0);

  cz( 4, 1) = -        sqrt( 5d0 /       43d0) *    5d0;
  cz( 4, 2) =          sqrt( 5d0 /       43d0) *   12d0;

  cz( 5, 1) =    2d0 * sqrt(15d0 /        7d0);
  cz( 6, 1) =    2d0 * sqrt(15d0 /        7d0);

  cz( 7, 1) = -  4d0 * sqrt(42d0 /     3685d0) *   14d0;
  cz( 7, 2) =    4d0 * sqrt(42d0 /     3685d0) *   25d0;
  cz( 8, 1) = -  4d0 * sqrt(42d0 /     3685d0) *   14d0;
  cz( 8, 2) =    4d0 * sqrt(42d0 /     3685d0) *   25d0;

  cz( 9, 1) =    4d0 * sqrt(10d0             ) /    3d0;

  cz(10, 1) =    4d0 * sqrt(70d0 /      103d0);

  cz(11, 1) =    3d0 / sqrt(        1072205d0) *  737d0;
  cz(11, 2) = -  3d0 / sqrt(        1072205d0) * 5140d0;
  cz(11, 3) =    3d0 / sqrt(        1072205d0) * 6020d0;

  cz(12, 1) = - 30d0 / sqrt(         492583d0) *  249d0;
  cz(12, 2) =   30d0 / sqrt(         492583d0) *  392d0;
  cz(13, 1) = - 30d0 / sqrt(         492583d0) *  249d0;
  cz(13, 2) =   30d0 / sqrt(         492583d0) *  392d0;

  cz(14, 1) =   10d0 * sqrt( 7d0 / 99258181d0) * 2970d0 / 3d0;
  cz(14, 2) = - 10d0 * sqrt( 7d0 / 99258181d0) * 5980d0 / 3d0;
  cz(14, 3) =   10d0 * sqrt( 7d0 / 99258181d0) * 5413d0 / 3d0;
  cz(15, 1) = - 10d0 * sqrt( 7d0 / 99258181d0) * 2970d0 / 3d0;
  cz(15, 2) =   10d0 * sqrt( 7d0 / 99258181d0) * 5980d0 / 3d0;
  cz(15, 3) =   10d0 * sqrt( 7d0 / 99258181d0) * 5413d0 / 3d0;

  cz(16, 1) =    2d0 * sqrt( 6d0 / 1089382547d0) *  70369d0;
  cz(16, 2) = -  2d0 * sqrt( 6d0 / 1089382547d0) * 322280d0;
  cz(16, 3) =    2d0 * sqrt( 6d0 / 1089382547d0) * 309540d0;
  cz(17, 1) =    2d0 * sqrt( 6d0 / 1089382547d0) *  70369d0;
  cz(17, 2) = -  2d0 * sqrt( 6d0 / 1089382547d0) * 322280d0;
  cz(17, 3) =    2d0 * sqrt( 6d0 / 1089382547d0) * 309540d0;

  cz(18, 1) = -  4d0 * sqrt(385d0 / 295894589d0) *   3322d0;
  cz(18, 2) =    4d0 * sqrt(385d0 / 295894589d0) *   4365d0;

  cz(19, 1) = -  4d0 * sqrt( 5d0 / 97d0) * 22d0;
  cz(19, 2) =    4d0 * sqrt( 5d0 / 97d0) * 35d0;

  cz(20, 1) = -  2.17600248d0;
  cz(20, 2) =   13.23551876d0;
  cz(20, 3) = - 16.15533716d0;
  cz(20, 4) =    5.95928883d0;

  cz(21, 1) =    2.17600248d0;
  cz(21, 2) = - 13.23551876d0;
  cz(21, 3) =   16.15533716d0;
  cz(21, 4) =    5.95928883d0;

  cz(22, 1) = -  2.47059083d0;
  cz(22, 2) =   33.14780774d0;
  cz(22, 3) = - 93.07966445d0;
  cz(22, 4) =   70.01749250d0;

  cr2  = cr .* cr;
  cr3  = cr .* cr2;
  cr4  = cr .* cr3;
  cr5  = cr .* cr4;
  cr6  = cr .* cr5;

  for iz = 1 : size(zc, 1)
    switch izp(iz)
      case  1,
        zer  = zer + zc(iz);
      case  2,
        zer  = zer + zc(iz) * (cz( 2,1) * cr) .* cos(ct);
      case  3,
        zer  = zer + zc(iz) * (cz( 3,1) * cr) .* sin(ct);
      case  4,
        zer  = zer + zc(iz) * (cz( 4,1) + cz(4,2) * cr2);
      case  5,
        zer  = zer + zc(iz) * (cz( 5,1) * cr2) .* sin(2 * ct);
      case  6,
        zer  = zer + zc(iz) * (cz( 6,1) * cr2) .* cos(2 * ct);
      case  7,
        zer  = zer + zc(iz) * (cz( 7,1) * cr + cz( 7,2) * cr3) .* sin(ct);
      case  8,
        zer  = zer + zc(iz) * (cz( 8,1) * cr + cz( 8,2) * cr3) .* cos(ct);
      case  9,
        zer  = zer + zc(iz) * (cz( 9,1) * cr3) .* sin(3 * ct);
      case 10,
        zer  = zer + zc(iz) * (cz(10,1) * cr3) .* cos(3 * ct);
      case 11,
        zer  = zer + zc(iz) ...
             * (cz(11,1) + cz(11,2) * cr2 + cz(11,3) * cr4);
      case 12,
        zer  = zer + zc(iz) ...
             * (cz(12,1) * cr2 + cz(12,2) * cr4) .* cos(2 * ct);
      case 13,
        zer  = zer + zc(iz) ...
             * (cz(13,1) * cr2 + cz(13,2) * cr4) .* sin(2 * ct);
      case 14,
        zer  = zer + zc(iz) ...
             * ((cz(14,1) * cr2 + cz(14,2) * cr4) .* cos(2 * ct) ...
             + cz(14,3) * cr4 .* cos(4 * ct));
      case 15,
        zer  = zer + zc(iz) ...
             * ((cz(15,1) * cr2 + cz(15,2) * cr4) .* sin(2 * ct) ...
             + cz(15,3) * cr4 .* sin(4 * ct));
      case 16,
        zer  = zer + zc(iz) ...
             * (cz(16,1) * cr + cz(16,2) * cr3 + cz(16,3) * cr5) .* cos(ct);
      case 17,
        zer  = zer + zc(iz) ...
             * (cz(17,1) * cr + cz(17,2) * cr3 + cz(17,3) * cr5) .* sin(ct);
      case 18,
        zer  = zer + zc(iz) ...
             * (cz(18,1) * cr3 + cz(18,2) * cr5) .* cos(3 * ct);
      case 19,
        zer  = zer + zc(iz) ...
             * (cz(19,1) * cr3 + cz(19,2) * cr5) .* sin(3 * ct);
      case 20,
        zer  = zer + zc(iz) ...
             * ((cz(20,1) * cr + cz(20,2) * cr3 + cz(20,3) * cr5) .* cos(ct) ...
             + cz(20,4) * cr5 .* cos(5 * ct));
      case 21,
        zer  = zer + zc(iz) ...
             * ((cz(21,1) * cr + cz(21,2) * cr3 + cz(21,3) * cr5) .* sin(ct) ...
             + cz(21,4) * cr5 .* sin(5 * ct));
      case 22,
        zer  = zer + zc(iz) ...
             * (cz(22,1) + cz(22,2) * cr2 + cz(22,3) * cr4 + cz(22,4) * cr6);
    end                 % switch izp(iz)
  end                   % for iz = 1 : size(zc, 1)

  abm  = zeros(ny, nx);
  abm(iay1 : iay2, iax1 : iax2) = zer;
end                     % function prop_hex_zernikes
