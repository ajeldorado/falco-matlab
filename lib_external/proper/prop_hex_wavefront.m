%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [bmo, ap, ph] = prop_hex_wavefront(bmi, nr, hr, hs, varargin)
%        [bmo, ap, ph] = prop_hex_wavefront(bmi, nr, hr, hs, varargin)
% Compute transmission and wavefront phase errors for a segmented
% hexagonal aperture and apply them to the current wavefront.  The
% hexagons making up the aperture have antialiased edges.  Each segment
% has user-defined hexagonal Zernike phase errors (up to Z22).
%
% Outputs:
% bmo  = beam structure output
% ap   = aperture transmission map
% ph   = phase map (microns of wavefront error)
%
% Required inputs:
% bmi  = beam structure input
% nr   = number of rings of hexagons in aperture
%        (e.g., 1 = a central hexagon surrounded by a ring of hexagons)
% hr   = distance from the center of a hexagonal segment to a vertex (m)
% hs   = distance between centers of adjacent hexagonal segments (m)
%
% Optional inputs:
% 'zernike_val'       = hexagonal Zernike polynomial coefficients
%                     = array of dimensions (nh, 22), where nh is the number
%                       of hexagonal segments (nh = nr * (nr + 1) * 3 + 1)
%                       Each line in this array specifies the coefficients
%                       for a segment.  The Zernikes are Noll-ordered
%                       (see prop_zernikes in the manual for a list of
%                       them or look at the prop_hex_zernikes.m routine
%                       in the PROPER library).
%                       The values are in meters of RMS wavefront error.
%                       Even if the center segment is made dark by
%                       'dark', there must be an entry for it.
% 'rotation'          = counter-clockwise rotation of the aperture about
%                       its center (degrees)
% 'darkcenter'        : if set, the central hexagonal segment
%                       transmission will be set to 0.0
% 'no_apply'          : if set, the aperture and phase maps are created
%                       but not applied to the wavefront
% 'xcenter'           = aperture center to wavefront center X (m)
% 'ycenter'           = aperture center to wavefront center Y (m)
%                       By default, the aperture is centered on the wavefront.

% 2007 Jul     jek  created idl routine
% 2016 Jul 26  gmg  Matlab translation
% 2017 Mar 16  gmg  Revised for keyword/value for optional inputs
% 2018 Jul 23  jek  Fixed segment numbering scheme (was not numbering central
%                   segment if it was set to dark by DARKCENTER)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  cx   = 0.0d0;         % aperture center to wavefront center X (m)
  cy   = 0.0d0;         % aperture center to wavefront center Y (m)
  dc   = 0;             % no dark center
  noap = 0;             % apply aperture to wavefront
  rot  = 0.0;           % counter-clockwise rotation of aperture (deg)
  zv   = 0.0;           % no hexagonal Zernike phase errors

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xcenter'}
        icav = icav + 1;
        cx   = varargin{icav};  % pixel coordinates of aperture center X
      case {'cy', 'ycenter'}
        icav = icav + 1;
        cy   = varargin{icav};  % pixel coordinates of aperture center Y
      case {'dc', 'darkcenter'}
        dc   = 1;               % dark center
      case {'noap', 'no_apply'}
        noap = 1;               % do not apply aperture to wavefront
      case {'rot', 'rotation'}
        icav = icav + 1;
        rot  = varargin{icav};  % counter-clockwise rotation of aperture (deg)
      case {'zv', 'zernike_val'}
        icav = icav + 1;
        zv   = varargin{icav};  % pixel coordinates of aperture center X
      otherwise
        error('prop_hex_wavefront: Unknown keyword: %s\n', varargin{icav});
    end
  end

  if zv == 0.0
    nh   = nr * (nr+1) * 3 + 1; % number of hexagonal segments
    zv   = zeros(nh, 22);       % Zernike coefficients (m)
  end
  ipe  = (sum(sum(abs(zv))) ~= 0);   % include phase errors
  zi   = [1 : 22];              % Zernike polynomial indices

  dx   = prop_get_sampling(bmi);% distance between pixels (m)
  [ny, nx] = size(bmi.wf);      % number of array points
  ap   = zeros(ny, nx);         % aperture array
  if ipe
    ph   = zeros(ny, nx);       % phase array (rad)
  else
    ph   = 0d0;                 % phase scalar
  end

  segi = 1;

  for ir = 0 : nr
    hx   = hs *  ir * cosd(30d0);       % hexagonal segment center X (m)
    hy   = hs * (ir * cosd(60d0) - nr); % hexagonal segment center Y (m)

    for iseg = 0 : 2 * nr - ir
% create hexagonal segment on one side
      hexx = cx + hx * cosd(rot) - hy * sind(rot);
      hexy = cy + hx * sind(rot) + hy * cosd(rot);
      seg  = prop_polygon(bmi, 6, hr, 'cx', hexx, 'cy', hexy, 'rot', rot);
      if ( (dc & (iseg == nr) & (ir == 0)) == 0 )
        ap   = ap + seg;
      end
% Remove apodization from seg
      seg  = (seg ~= 0);
      if ipe
          ph = ph + seg .* prop_hex_zernikes(zi, zv(segi,:)', nx, dx, hr, ...
                           'hx', hexx, 'hy', hexy, 'ang', rot);
      end
      segi = segi + 1;

% create hexagonal segment on opposite side
      if (ir ~= 0)
        hexx = cx - hx * cosd(rot) - hy * sind(rot);
        hexy = cy - hx * sind(rot) + hy * cosd(rot);
        seg  = prop_polygon(bmi, 6, hr, 'cx', hexx, 'cy', hexy, 'rot', rot);
        ap   = ap + seg;
% Remove apodization from seg
        seg  = (seg ~= 0);
        if ipe
          ph = ph + seg .* prop_hex_zernikes(zi, zv(segi,:)', nx, dx, hr, ...
                           'hx', hexx, 'hy', hexy, 'ang', rot);
        end
        segi = segi + 1;
      end

      hy   = hy + hs;
    end                 % for iseg = 0 : 2 * nr - ir
  end                   % for ir = 0 : nr

  if ~noap
    bmo  = prop_multiply(bmi, ap);
    if ipe
      bmo = prop_add_phase(bmo, ph);
    end
  else
    bmo  = bmi;
  end
end                     % prop_hex_wavefront
