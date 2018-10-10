%   Copyright 2016, 2017, 2018 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt
%   Centering option by A.J. Riggs
% -------------------------------------------------------------------------


function [bmo, amap] = propcustom_zernikes(bmi, zn, zv, varargin)
%        [bmo, amap] = prop_zernikes(bmi, zn, zv, varargin)
% Add Zernike-polynomial wavefront errors to current wavefront.
% Noll ordering is used and a circular system is assumed.
% An arbitrary number of Zernikes normalized for an unobscured circular
% region can be computed, but only the first 22 Zernikes can be computed
% normalized for a centrally-obscured region.
% The user specifies 1D arrays containing the Zernike polynomial
% coefficient indices, the respective coefficients, and if an obstructed
% circular aperture the central obscuration ratio.  A wavefront error
% map will be generated and added to the current wavefront.
%
% Outputs:
% bmo  = beam wavefront structure output
% amap = aberration map that was added to the wavefront (in meters rms)
%
% Required inputs:
% bmi  = beam wavefront structure input
% zn   = 1D array specifying which Zernike polynomials to include
% zv   = 1D array containing Zernike coefficients (in meters of RMS
%        wavefront phase error or dimensionless RMS amplitude error)
%        for Zernike polynomials indexed by "zn"
%
% Optional inputs:
% 'eps'               = central obscuration ratio (0.0-1.0); default 0.0
% 'amplitude'         : specifies that the Zernike values in "zv"
%                       represent wavefront RMS amplitude (rather than
%                       phase) variation.  The current wavefront will be
%                       multiplied by the generated map.
% 'name'              = string containing name of surface that will be
%                       printed when executed
% 'no_apply'          : If set, the aberration map will be generated,
%                       but will not be applied to the wavefront.  This
%                       is useful if you just want to generate a map for
%                       your own use and modification (e.g., to create
%                       an error map for a multi-segmented system, each
%                       with its own aberration map).
% 'radius'            = The radius to which the Zernike polynomials are
%                       normalized.  If not specified, the pilot beam
%                       radius is used. (m)
%
% Zernike index and corresponding aberration for 1st 22 zernikes
%        1 : Piston
%        2 : X tilt
%        3 : Y tilt
%        4 : Focus
%        5 : 45 degree astigmatism
%        6 : 0 degree astigmatism
%        7 : Y coma
%        8 : X coma
%        9 : Y clover (trefoil)
%       10 : X clover (trefoil)
%       11 : 3rd order spherical
%       12 : 5th order 0 degree astig
%       13 : 5th order 45 degree astig
%       14 : X quadrafoil
%       15 : Y quadrafoil
%       16 : 5th order X coma
%       17 : 5th order Y coma
%       18 : 5th order X clover
%       19 : 5th order Y clover
%       20 : X pentafoil
%       21 : Y pentafoil
%       22 : 5th order spherical

% 2005 Feb     jek  created idl routine
% 2016 Jun 23  gmg  Matlab translation
% 2017 Mar 10  gmg  Revised for keyword/value for optional inputs
% 2018 Aug 10  a_r  Added bmi.centering option for interpixel centering.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  ampl = 0;             % Zernike coefficients represent RMS phase variation
  bc   = 0.0;           % central obscuration ratio
  brad = prop_get_beamradius(bmi);
  noap = 0;             % apply aberration map
  sfnm = '';            % surface name

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'ampl', 'amplitude'}
        ampl = 1;               % Zernike coefficients represent RMS amplitude
      case {'bc', 'eps'}
        icav = icav + 1;
        bc   = varargin{icav};  % input array center X (pixels)
      case {'zr', 'radius'}
        icav = icav + 1;
        brad = varargin{icav};  % radius to normalize Zernike polynomials
      case {'noap', 'no_apply'}
        noap = 1;               % don't apply aberration error map
      case {'sfnm', 'name'}
        icav = icav + 1;
        sfnm = varargin{icav};  % surface name
      otherwise
        error('prop_zernikes: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  propcommon;

  if print_it & ~noap
    if strcmp(sfnm, '');
      fprintf(1, 'Applying aberrations\n');
    else
      fprintf(1, 'Applying aberrations at %s\n', sfnm);
    end
  end

  maxz = max(zn);               % maximum Zernike index

  if (bc ~= 0.0) & (maxz > 22)
    error('PROP_ZERNIKES:', ...
      'Maximum index for an obscured Zernike polynomial is 22.\n');
  end

  [ny, nx] = size(bmi.wf);
  amap = zeros(ny, nx);         % aberration map added to wavefront (m)
  
  
%   cx   = fix(nx / 2) + 1;       % center pixel x (pixels)
%   cy   = fix(ny / 2) + 1;       % center pixel y (pixels)
%   [px, py] = meshgrid((1 : nx) - cx, (1 : ny) - cy); % Create coordinate arrays (pixels)

% Create coordinate arrays, pixel-centered (units of pixels)
xs = -nx/2:(nx/2-1);
ys = -ny/2:(ny/2-1);
[px, py] = meshgrid(xs,ys);

% Create coordinate arrays, interpixel-centered (units of pixels)
  if(isfield(bmi,'centering'))
      if(strcmpi(bmi.centering,'interpixel'))
          fprintf('Shifting Zernikes to have interpixel centering.\n')
            xs = -(nx-1)/2:(nx-1)/2;
            ys = -(ny-1)/2:(ny-1)/2;
            [px, py] = meshgrid(xs,ys);
      end
  end
  
  r    = sqrt(px.^2 + py.^2) * bmi.dx / brad;   % normalized radius
  t    = atan2(py, px);                         % azimuth angle (radians)

  if bc == 0.0                  % Zernikes for unobscured region

% Get list of executable equations defining Zernike polynomials
    zca  = prop_noll_zernikes(maxz);

    for iz = 1 : max(size(zn))
      ab   = eval(char(zca(zn(iz))));
      amap = amap + zv(iz) * ab;
    end                         % for iz = 1 : max(size(zn))

  else                          % Zernikes for a centrally-obscured region
    sr02 = sqrt( 2);
    sr03 = sqrt( 3);
    sr05 = sqrt( 5);
    sr06 = sqrt( 6);
    sr07 = sqrt( 7);
    sr10 = sqrt(10);

    r2   = r .* r;
    r3   = r .* r2;
    r4   = r .* r3;
    r5   = r .* r4;
    r6   = r .* r5;

    for iz = 1 : max(size(zn))
      switch zn(iz)
        case  1
          ab = 1;
        case  2
          ab = 2        * cos(  t) .* r  ...
                            / sqrt(bc^2 + 1);
        case  3
          ab = 2        * sin(  t) .* r  ...
                            / sqrt(bc^2 + 1);
        case  4
          ab =     sr03 * (1 + bc^2 - 2*r2) ...
             / (bc^2 - 1);
        case  5
          ab =     sr06 * sin(2*t) .* r2 ...
                            / sqrt(bc^4 + bc^2 + 1);
        case  6
          ab =     sr06 * cos(2*t) .* r2 ...
                            / sqrt(bc^4 + bc^2 + 1);
        case  7
          ab = 2 * sr02 * sin(  t) .* r  ...
             .* (2 + 2*bc^4 - 3*r2 + (2 - 3*r2) * bc^2) ...
             / (bc^2 - 1)   / sqrt(bc^6 + 5*bc^4 + 5*bc^2 + 1);
        case  8
          ab = 2 * sr02 * cos(  t) .* r  ...
             .* (2 + 2*bc^4 - 3*r2 + (2 - 3*r2) * bc^2) ...
             / (bc^2 - 1)   / sqrt(bc^6 + 5*bc^4 + 5*bc^2 + 1);
        case  9
          ab = 2 * sr02 * sin(3*t) .* r3 ...
                            / sqrt(bc^6 + bc^4 + bc^2 + 1);
        case 10
          ab = 2 * sr02 * cos(3*t) .* r3 ...
                            / sqrt(bc^6 + bc^4 + bc^2 + 1);
        case 11
          ab =     sr05 * (1 + bc^4 - 6*r2 + 6*r4 + (4 - 6*r2)*bc^2) ...
             / (bc^2 - 1)^2;
        case 12
          ab =     sr10 * cos(2*t) .* r2 .* (3 + 3*bc^6 - 4*r2 ...
             + (3 - 4*r2) * (bc^4 + bc^2)) ...
             / (bc^2 - 1)   / sqrt(bc^4 + bc^2 + 1) ...
             / sqrt(bc^8 + 4*bc^6 + 10*bc^4 + 4*bc^2 + 1);
        case 13
          ab =     sr10 * sin(2*t) .* r2 .* (3 + 3*bc^6 - 4*r2 ...
             + (3 - 4*r2) * (bc^4 + bc^2)) ...
             / (bc^2 - 1)   / sqrt(bc^4 + bc^2 + 1) ...
             / sqrt(bc^8 + 4*bc^6 + 10*bc^4 + 4*bc^2 + 1);
        case 14
          ab =     sr10 * cos(4*t) .* r4 ...
                            / sqrt(bc^8 + bc^6 + bc^4 + bc^2 + 1);
        case 15
          ab =     sr10 * sin(4*t) .* r4 ...
                            / sqrt(bc^8 + bc^6 + bc^4 + bc^2 + 1);
        case 16
          ab = 2 * sr03 * cos(  t) .* r  .* (3 + 3*bc^8 - 12*r2 ...
             + 10*r4 - 12*(r2 - 1)*bc^6 + 2*(15 - 24*r2 + 5*r4)*bc^4 ...
             + 4*(3 - 12*r2 + 10*r4)*bc^2) ...
             / (bc^2 - 1)^2 / sqrt(bc^4 + 4*bc^2 + 1) ...
             / sqrt(bc^6 + 9*bc^4 + 9*bc^2 + 1);
        case 17
          ab = 2 * sr03 * sin(  t) .* r  .* (3 + 3*bc^8 - 12*r2 ...
             + 10*r4 - 12*(r2 - 1)*bc^6 + 2*(15 - 24*r2 + 5*r4)*bc^4 ...
             + 4*(3 - 12*r2 + 10*r4)*bc^2) ...
             / (bc^2 - 1)^2 / sqrt(bc^4 + 4*bc^2 + 1) ...
             / sqrt(bc^6 + 9*bc^4 + 9*bc^2 + 1);
        case 18
          ab = 2 * sr03 * cos(3*t) .* r3 .* (4 + 4*bc^8 - 5*r2 ...
             + (4 - 5*r2) * (bc^6 + bc^4 + bc^2))...
             / (bc^2 - 1)   / sqrt(bc^6 + bc^4 + bc^2 + 1) ...
             / sqrt(bc^12+4*bc^10+10*bc^8+20*bc^6+10*bc^4+4*bc^2+1);
        case 19
          ab = 2 * sr03 * sin(3*t) .* r3 .* (4 + 4*bc^8 - 5*r2 ...
             + (4 - 5*r2) * (bc^6 + bc^4 + bc^2))...
             / (bc^2 - 1)   / sqrt(bc^6 + bc^4 + bc^2 + 1) ...
             / sqrt(bc^12+4*bc^10+10*bc^8+20*bc^6+10*bc^4+4*bc^2+1);
        case 20
          ab = 2 * sr03 * cos(5*t) .* r5 ...
                            / sqrt(bc^10 + bc^8 + bc^6 + bc^4 + bc^2 + 1);
        case 21
          ab = 2 * sr03 * sin(5*t) .* r5 ...
                            / sqrt(bc^10 + bc^8 + bc^6 + bc^4 + bc^2 + 1);
        case 22
          ab =     sr07 * (1 + bc^6 - 12*r2 + 30*r4 - 20*r6 ...
             + (9 - 36*r2 + 30*r4) * bc^2 + (9 - 12*r2) * bc^4) ...
             / (bc^2 - 1)^3;
      end                       % switch zn(iz)
      amap = amap + zv(iz) * ab;
    end                         % for iz = 1 : max(size(zn))
  end

  bmo = bmi;
  if ~noap
    if ampl
      bmo.wf = bmi.wf .* prop_shift_center(amap);
    else
      bmo.wf = bmi.wf .* exp(1i * 2.0d0 * pi * prop_shift_center(amap) / bmi.wl);
    end
  end
end                     % function prop_zernikes
