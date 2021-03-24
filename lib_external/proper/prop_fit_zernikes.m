%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [zc, fit] = prop_fit_zernikes(wf, msk, apr, nz, varargin)
%        [zc, fit] = prop_fit_zernikes(wf, msk, apr, nz, varargin)
% Fit circular Zernike polynomials to a 2D error map.  The user provides
% the error map, a 2D mask (zero or one) that defines valid data values
% in the map, and the radius of the aperture.
%
% Outputs:
% zc  = Fitted Zernike coefficients, ordered such that zc(1) is the
%       first Zernike polynomial (Z1, piston).
%       The units are the same as those in wf.
% fit = 2D map of the Zernike polynomials fitted to wf.
%       This map can be directly subtracted from wf, for instance.
%
% Required inputs:
% wf  = A 2D array containing the aberration map.  The returned Zernike
%       coefficients will have the same data units as this map.
% msk = 2D mask array indicating which corresponding values in wf are
%       valid.  A value of 1 indicates a valid point, 0 is a bad point
%       or is obscured.
% apr = Aperture radius of the beam in the map in pixels.  The Zernike
%       polynomials are normalized for a circle of this radius.
% nz  = Maximum number of Zernike polynomials to fit (1 to nz using the
%       Noll ordering scheme).  This is arbitrary for unobscured
%       polynomials.  For obscured ones, the max allowed is 22.
%
% Optional inputs:
% 'obscuration_ratio' = Obscuration ratio of the central obscuration
%                       radius to the aperture radius.  Specifying this
%                       value will cause Zernike polynomials normalized
%                       for an obscured aperture to be fit rather than
%                       unobscured Zernikes, which is the default.
% 'xc'                = Specifies the center of the wavefront in the
% 'yc'                  wavefront array in pixels, with the center of
%                       the lower-left pixel being (0.0, 0.0).  By
%                       default, the wavefront center is at the center
%                       of the array (fix(nx/2) + 1, fix(ny/2) + 1).

% 2005 Feb     jek  created idl routine
% 2016 Mar 28  gmg  Matlab translation
% 2017 Jan 03  gmg  Switched to explicit least squares
%                   Fixed bug in which fitted Zernikes map was returned
%                   at shrunken wavefront sampling
% 2017 Jan 25  gmg  Removed shr option; output "fit" map is not masked
% 2017 Feb 10  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bc   = 0.0;                   % obscuration ratio
  [ny, nx] = size(wf);
  cx = fix(nx / 2) + 1;         % center pixel x (pixels)
  cy = fix(ny / 2) + 1;         % center pixel y (pixels)

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'eps', 'bc', 'obscuration_ratio'}
        icav = icav + 1;
        bc   = varargin{icav};
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};
      otherwise
        error('prop_fit_zernikes: Unknown keyword: %s\n', varargin{icav});
    end
  end

  if bc > 0.0 & nz > 22
    error('Proper:PROP_FIT_ZERNIKES', ...
          'Limited to first 22 Obscured Zernikes.\n');
  end

  for ik = 0 : nargout - 1
% When ik = 0, fit the "good" values in the wavefront.
% If the "fit" output parameter is requested, go through again (ik = 1),
% this time creating the aberration map for all pixels using
% the fitted coefficients.

% Create coordinate arrays (pixels)
    [px, py] = meshgrid((1 : nx) - cx, (1 : ny) - cy);
    r    = sqrt(px.^2 + py.^2) / apr;   % normalized radius
    t    = atan2(py, px);               % azimuth angle (radians)

% During fitting stage, only use pixels in pupil region
    if ik == 0
      imr  = (msk ~= 0) & (r <= 1.0);
      r    = r(imr);            % column array of normalized radius
      t    = t(imr);            % column array of azimuth angle (radians)
      wf   = wf(imr);           % column array of wavefront
    else
      r    = r(:);              % column array of normalized radius
      t    = t(:);              % column array of azimuth angle (radians)
    end                         % if ik == 0
    nrca = size(r, 1);          % number of rows in column arrays

    ab   = zeros(nrca, nz);     % 1D array values for each Zernike

    if bc > 0.0
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
        ab(:, 1) = 1.0;
      if nz >  1
        ab(:, 2) = 2        * cos(  t) .* r  ...
                                / sqrt(bc^2 + 1);
      end
      if nz >  2
        ab(:, 3) = 2        * sin(  t) .* r  ...
                                / sqrt(bc^2 + 1);
      end
      if nz >  3
        ab(:, 4) =     sr03 * (1 + bc^2 - 2*r2) ...
                 / (bc^2 - 1);
      end
      if nz >  4
        ab(:, 5) =     sr06 * sin(2*t) .* r2 ...
                                / sqrt(bc^4 + bc^2 + 1);
      end
      if nz >  5
        ab(:, 6) =     sr06 * cos(2*t) .* r2 ...
                                / sqrt(bc^4 + bc^2 + 1);
      end
      if nz >  6
        ab(:, 7) = 2 * sr02 * sin(  t) .* r  ...
                 .* (2 + 2*bc^4 - 3*r2 + (2 - 3*r2) * bc^2) ...
                 / (bc^2 - 1)   / sqrt(bc^6 + 5*bc^4 + 5*bc^2 + 1);
      end
      if nz >  7
        ab(:, 8) = 2 * sr02 * cos(  t) .* r  ...
                 .* (2 + 2*bc^4 - 3*r2 + (2 - 3*r2) * bc^2) ...
                 / (bc^2 - 1)   / sqrt(bc^6 + 5*bc^4 + 5*bc^2 + 1);
      end
      if nz >  8
        ab(:, 9) = 2 * sr02 * sin(3*t) .* r3 ...
                                / sqrt(bc^6 + bc^4 + bc^2 + 1);
      end
      if nz >  9
        ab(:,10) = 2 * sr02 * cos(3*t) .* r3 ...
                                / sqrt(bc^6 + bc^4 + bc^2 + 1);
      end
      if nz > 10
        ab(:,11) =     sr05 * (1 + bc^4 - 6*r2 + 6*r4 + (4 - 6*r2)*bc^2) ...
                 / (bc^2 - 1)^2;
      end
      if nz > 11
        ab(:,12) =     sr10 * cos(2*t) .* r2 .* (3 + 3*bc^6 - 4*r2 ...
                 + (3 - 4*r2) * (bc^4 + bc^2)) ...
                 / (bc^2 - 1)   / sqrt(bc^4 + bc^2 + 1) ...
                 / sqrt(bc^8 + 4*bc^6 + 10*bc^4 + 4*bc^2 + 1);
      end
      if nz > 12
        ab(:,13) =     sr10 * sin(2*t) .* r2 .* (3 + 3*bc^6 - 4*r2 ...
                 + (3 - 4*r2) * (bc^4 + bc^2)) ...
                 / (bc^2 - 1)   / sqrt(bc^4 + bc^2 + 1) ...
                 / sqrt(bc^8 + 4*bc^6 + 10*bc^4 + 4*bc^2 + 1);
      end
      if nz > 13
        ab(:,14) =     sr10 * cos(4*t) .* r4 ...
                                / sqrt(bc^8 + bc^6 + bc^4 + bc^2 + 1);
      end
      if nz > 14
        ab(:,15) =     sr10 * sin(4*t) .* r4 ...
                                / sqrt(bc^8 + bc^6 + bc^4 + bc^2 + 1);
      end
      if nz > 15
        ab(:,16) = 2 * sr03 * cos(  t) .* r  .* (3 + 3*bc^8 - 12*r2 ...
                 + 10*r4 - 12*(r2 - 1)*bc^6 + 2*(15 - 24*r2 + 5*r4)*bc^4 ...
                 + 4*(3 - 12*r2 + 10*r4)*bc^2) ...
                 / (bc^2 - 1)^2 / sqrt(bc^4 + 4*bc^2 + 1) ...
                 / sqrt(bc^6 + 9*bc^4 + 9*bc^2 + 1);
      end
      if nz > 16
        ab(:,17) = 2 * sr03 * sin(  t) .* r  .* (3 + 3*bc^8 - 12*r2 ...
                 + 10*r4 - 12*(r2 - 1)*bc^6 + 2*(15 - 24*r2 + 5*r4)*bc^4 ...
                 + 4*(3 - 12*r2 + 10*r4)*bc^2) ...
                 / (bc^2 - 1)^2 / sqrt(bc^4 + 4*bc^2 + 1) ...
                 / sqrt(bc^6 + 9*bc^4 + 9*bc^2 + 1);
      end
      if nz > 17
        ab(:,18) = 2 * sr03 * cos(3*t) .* r3 .* (4 + 4*bc^8 - 5*r2 ...
                 + (4 - 5*r2) * (bc^6 + bc^4  + bc^2))...
                 / (bc^2 - 1)   / sqrt(bc^6 + bc^4 + bc^2 + 1) ...
                 / sqrt(bc^12+4*bc^10+10*bc^8+20*bc^6+10*bc^4+4*bc^2+1);
      end
      if nz > 18
        ab(:,19) = 2 * sr03 * sin(3*t) .* r3 .* (4 + 4*bc^8 - 5*r2 ...
                 + (4 - 5*r2) * (bc^6 + bc^4  + bc^2))...
                 / (bc^2 - 1)   / sqrt(bc^6 + bc^4 + bc^2 + 1) ...
                 / sqrt(bc^12+4*bc^10+10*bc^8+20*bc^6+10*bc^4+4*bc^2+1);
      end
      if nz > 19
        ab(:,20) = 2 * sr03 * cos(5*t) .* r5 ...
                                / sqrt(bc^10 + bc^8 + bc^6 + bc^4 + bc^2 + 1);
      end
      if nz > 20
        ab(:,21) = 2 * sr03 * sin(5*t) .* r5 ...
                                / sqrt(bc^10 + bc^8 + bc^6 + bc^4 + bc^2 + 1);
      end
      if nz > 21
        ab(:,22) =     sr07 * (1 + bc^6 - 12*r2 + 30*r4 - 20*r6 ...
                 + (9 - 36*r2 + 30*r4) * bc^2 + (9 - 12*r2) * bc^4) ...
                 / (bc^2 - 1)^3;
      end
    else
      zca  = prop_noll_zernikes(nz);    % cell array of Zernike equations
      for iz = 1 : nz
        ab(:,iz) = eval(char(zca(iz)));
      end
    end                 % if bc > 0.0

    if ik == 0          % fit Zernikes to wavefront
      zc   = ab \ wf;   % least squares fit
    else
% Create map of fitted Zernikes at non-shrunken scale
      fit  = zeros(nrca, 1);
      for iz = 1 : nz
        fit  = fit  + zc(iz) * ab(:, iz);
      end
      fit  = reshape(fit, [ny, nx]);
    end                 % if ik == 0
  end                   % for ik = 0 : nargout - 1
end                     % function prop_fit_zernikes
