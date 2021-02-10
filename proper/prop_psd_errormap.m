%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [bm, map] = prop_psd_errormap(bm, resf, clp, fple, varargin)
%        [bm, map] = prop_psd_errormap(bm, resf, clp, fple, varargin)
% Create a realization of a two-dimensional surface, wavefront, or
% amplitude error map for a specified power-spectral-density (PSD)
% profile.  This map is applied to the current wavefront.
%
% Outputs:
% bm   = output wavefront structure
% map  = wavefront (meters), surface (meters), or amplitude map
%
% Required inputs:
% bm   = the current wavefront structure
% resf = low spatial frequency RMS error per spatial frequency
%        By default, this is the wavefront error (meters^4).
%        If the "mirror" switch is set, this is assumed to be the RMS
%        surface error (meters^4).  If the "rms" switch is set, then the
%        entire error map will be renormalized to have an RMS of "resf"
%        (if the PSD is vastly dominated by very low spatial frequency
%        errors, then the map in this case may not have the specified
%        RMS over the area of the beam but will over the grid, leading
%        to unexpected errors).
%        If "amplitude" is specified, then "resf" is assumed to be the RMS
%        amplitude error of the entire map (the rms switch is ignored).
% clp  = correlation length parameter (cycles/meter)
%        This basically indicates where the PSD curve transitions from a
%        horizontal line at low spatial frequencies to a sloping line at
%        high ones.
% fple = high frequency falloff power law exponent
%
% Optional inputs:
% 'file'              = the filename containing the error map generated
%                       by this routine.
%       This is used in one of two ways:
%   1) If the file exists, it is read in and that map is used instead of
%      creating a new map.  In this case, the parameters on the command
%      line are ignored, except when "ampl" is specified, in which case
%      the map read from the file will be assumed to be an amplitude
%      error map and will be adjusted to have a mean value specified by
%      the value of "ampl" over the entire array.
%      There are no checks made to verify the the command line parameters
%      and those used to generate the error map in the file are the same.
%   2) If the file doesn't exist, a map is generated and written to that
%      file.
%
% 'amplitude'         = an amplitude error map will be created with a
%                       maximum amplitude of "amplitude".  Remember that
%                       intensity is the square of the amplitude, so
%                       specifying ampl = 0.90 will result in a maximum
%                       intensity transmission of 0.81.  Because the map
%                       is normalized by the maximum value in the entire
%                       array, the maximum amplitude may be lower within
%                       the illuminated region of the beam.
%
% 'mirror'            : indicates that the specified PSD is for the
%                       surface error on a mirror, so the resulting
%                       wavefront error will be twice the map error.
%                       The map returned in "map" will be surface, not
%                       wavefront, error.
%
% 'no_apply'          : indicates a map will be generated but not
%                       applied to the wavefront (added if wavewfront or
%                       surface, multiplied if amplitude).
%                       This is useful if you wish to use the map for
%                       some custom purpose.
%
% 'rms'               : indicates that the error map array is to be
%                       normalized to have an RMS value about the mean
%                       specified by "resf".  Ignored if "amplitude" is
%                       set.  By default, "resf" specifies the RMS
%                       amplitude of just the low-spatial-frequency
%                       error component of the PSD.
%                       WARNING: The map is renormalized only during the
%                       initial creation - if a map of the same name was
%                       previously created and so is read in from a file,
%                       the map in the file will not be renormalized.
%
% 'inclination'       = inclination in the Y-axis (degrees)
% 'rotation'          = rotation about the Z-axis (degrees)
%                       Specify the angles of the surface or wavefront
%                       error plane relative to the incoming beam.  This
%                       approximately accounts for the projection of the
%                       beam onto the inclined surface and the resulting
%                       difference in spatial frequency scales of errors
%                       along the wavefront axes.  See the PSD section
%                       of the manual for more info.
%
% 'tpf'               : indicates that the Terrestrial Planet Finder
%                       2D PSD shape is to be used:
%
%                     resf
%   PSD_2D(k) = ----------------
%                           fple
%                1 + (k/clp)
%
%    where k and clp are in cycles/meter.
%    The default PSD shape is used if "tpf" is not specified:
%
%                         resf
%   PSD_2D(k) = ------------------------
%                             (fple+1)/2
%               (1 + (k/clp)^2)
%
%    This is the K-correlation form.
%    (see Church et al., Proc. of the SPIE, v. 1165, 136 (1989)).
%
% 'max_frequency'     = maximum spatial frequency (cycles/meter) in the
%                       generated map.  This can be used to prevent high
%                       spatial frequency components from generating
%                       aliasing errors when a map is resampled.

% 2005 Feb     jek  created idl routine
% 2007 Jun     jek  changed the behavior of the AMPLITUDE keyword; it now
%                   specifies the maximum amplitude (sqrt of intensity)
%                   of the error map rather than the mean amplitude.
% 2014 Aug 26  gmg  Matlab translation
% 2016 Dec 20  gmg  Added maxf option
% 2017 Apr 10  gmg  Revised for keyword/value for optional inputs
% 2017 Nov 01  gmg  Fixed bug: if 'file' isn't given, don't create map file
%                   Changed rng to use different seed each time
%                   prop_psd_errormap.m is run
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  global print_it
% Set default values of input parameters
  ampl = 0.0;           % not an amplitude error map
  angy = 0.0;           % inclination in the Y-axis (degrees)
  flnm = '';            % no pre-existing map file
  maxf = 0.0;           % maximum spatial frequency (cycles / meter)
  mirr = 0;             % error map is not surface error on a mirror
  noap = 0;             % apply map to wavefront
  rms  = 0;             % do not normalize map rms about mean = resp
  rotz = 0.0;           % rotation about the Z-axis (degrees)
  tpf  = 0;             % use default PSD shape

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'ampl', 'amplitude'}
        icav = icav + 1;
        ampl = varargin{icav};  % maximum amplitude of amplitude error map
      case {'angy', 'inclination'}
        icav = icav + 1;
        angy = varargin{icav};  % inclination in the Y-axis (deg)
      case {'file'}
        icav = icav + 1;
        flnm = varargin{icav};  % pre-existing map file name
      case {'maxf', 'max_frequency'}
        icav = icav + 1;
        maxf = varargin{icav};  % maximum spatial frequency
      case {'mirr', 'mirror'}
        mirr = 1;               % error map is surface error on a mirror
      case {'noap', 'no_apply'}
        noap = 1;               % do not apply map to wavefront
      case {'rms'}
        rms  = 1;               % normalize map rms about mean = resp
      case {'rotz', 'rotation'}
        icav = icav + 1;
        rotz = varargin{icav};  % rotation about the Z-axis (deg)
      case {'tpf'}
        tpf  = 1;               % use Terrestrial Planet Finder 2D PSD shape
      otherwise
        error('prop_psd_errormap: Unknown keyword: %s\n', varargin{icav});
    end
  end

  if exist(flnm) ~= 2                   % file does not exist, create map

    if print_it
      if ampl == 0.0                    % create phase aberration map
        fprintf(1, '  Creating phase aberration map from PSD\n');
      else                              % create ampli aberration map
        fprintf(1, '  Creating amplitude aberration map from PSD\n');
      end
    end

    [ny, nx] = size(bm.wf);
    dx   = prop_get_sampling(bm);       % array sample spacing x (m)
    rx   = nx * dx;                     % array range x (m)
    dk   = 1.0 / rx;                    % cycles / meter

    icpx = fix(nx / 2) + 1;             % index of center of array x
    icpy = fix(ny / 2) + 1;             % index of center of array y
    kx   = [1 : nx] - icpx;
    ky   = [1 : ny] - icpy;
    [ax, ay] = meshgrid(kx, ky);
    ax   = ax * cosd(-rotz) - ay * sind(-rotz);
    ay   = ax * sind(-rotz) + ay * cosd(-rotz);
    ay   = ay * cosd(-angy);
    kpsd = dk * sqrt(ax.^2 + ay.^2);    % cycles / meter

    if tpf == 1                         % use TPF 2D PSD shape
      psd2 = resf ./ (1.0 + (kpsd / clp)     .^   fple);
    else
      psd2 = resf ./ (1.0 + (kpsd / clp).^2) .^ ((fple + 1.0) / 2.0);
    end
    psd2(icpy, icpx) = 0.0;             % no piston
    rpsd = dk * sqrt(sum(sum(psd2)));   % RMS error from PSD

% Create realization of PSD using random phases
    rng('shuffle');             % set Mersenne Twistor generator state
%   fixe = zeros(ny, nx);
%   fixe(icpy, icpx) = 1.0;
%   cphs = exp(complex(0.0, 2.0 * pi * fixe(ny, nx) - pi));
%   cphs = exp(complex(0.0, pi * ones(ny, nx)));
    cphs = exp(complex(0.0, 2.0 * pi * rand(ny, nx) - pi));

% If specified, zero-out spatial frequencies > maxf
    if maxf > 0.0
      psd2 = psd2 .* (kpsd <= maxf);
    end

    map  = real(fft2(prop_shift_center(cphs .* sqrt(psd2) / dk, 1))) / rx^2;
% Note that Matlab's fft has the 1/N normalization in the inverse fft

% Force realized map to have RMS expected from PSD
    rmap = sqrt(sum(sum(map.^2)) / nx / ny);

    if (ampl ~= 0.0) | (rms == 1)
      map  = map * resf / rmap;
    else
      map  = map * rpsd / rmap;
    end

  else                          % map exists; read it in

    if print_it
      if ampl == 0.0                    % read in phase aberration map
        fprintf(1, '  PSD-realization phase map exists.');
        fprintf(1, '  Reading in %s\n', flnm);
      else                              % read in ampli aberration map
        fprintf(1, '  PSD-realization amplitude map exists.');
        fprintf(1, '  Reading in %s\n', flnm);
      end
    end

    map  = prop_readmap(bm, flnm);

  end                           % if exist(flnm) ~= 2

  if ampl ~= 0.0
    type = 'amplitude';
    mmap = max(max(real(map)));
    map  = map + ampl - mmap;
    if noap ~= 1
      bm.wf = bm.wf .* map;
    end
  elseif mirr == 1
    type = 'mirror surface';
    if noap ~= 1
      bm.wf = bm.wf .* exp(complex(0.0, 4.0 * pi * map / bm.wl));
    end
  else
    type = 'wavefront';
    if noap ~= 1
      bm.wf = bm.wf .* exp(complex(0.0, 2.0 * pi * map / bm.wl));
    end
  end

  map  = prop_shift_center(map);

  if (~strcmp(flnm,'')&exist(flnm)~=2)  % file does not exist, create it
    fprintf(1, '  Writing out PSD realization map to %s\n', flnm);

    import matlab.io.*;
    fptr = fits.createFile(flnm);
    fits.createImg(fptr, 'double_img', [ny nx]);
    fits.deleteKey(fptr, 'EXTEND');
    fits.writeDate(fptr);
    fits.writeImg(fptr, map);

    fits.writeKey(fptr, 'MAPTYPE', type, ' error map type');
    if tpf == 1
      fits.writeKey(fptr, 'PSDTYPE', 'TPF');
    end
    fits.writeKey(fptr, 'X_UNIT', 'meters', ' X-Y units');
    fits.writeKey(fptr, 'PIXSIZE', dx, ' spacing in meters');
    if ~strcmp(type, 'amplitude')
      fits.writeKey(fptr, 'Z_UNIT', 'meters', ' Error units');
      fits.writeKey(fptr, 'PSD_AMP', resf, ...
        ' PSD low frequency RMS amplitude (m^4)');
    else
      fits.writeKey(fptr, 'PSD_AMP', resf, ...
        ' PSD low frequency RMS amplitude (amp^2 m^2)');
    end
    fits.writeKey(fptr, 'PSD_B', clp, ' PSD correlation length (cycles/m)');
    fits.writeKey(fptr, 'PSD_C', fple, ' PSD high frequency power law');
    icpx = int16(fix(nx / 2) + 1);      % index of center of array x
    icpy = int16(fix(ny / 2) + 1);      % index of center of array y
    fits.writeKey(fptr, 'XC_PIX', icpx, ' Center X pixel coordinate');
    fits.writeKey(fptr, 'YC_PIX', icpy, ' Center Y pixel coordinate');
    if maxf > 0.0
      fits.writeKey(fptr, 'MAXFREQ', maxf, ...
        ' Maximum spatial frequency in cycles/meter');
    end
    fits.closeFile(fptr);
  end                   % if (~strcmp(flnm,'')&exist(flnm)~=2)
end                     % function prop_psd_errormap
