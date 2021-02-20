%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [bm, map] = prop_errormap(bm, flnm, varargin)
%        [bm, map] = prop_errormap(bm, flnm, varargin)
% Read in an amplitude, surface, or wavefront error map from a
% Flexible Image Transport System (FITS) file.
% One of the types (amplitude, mirror surface, or wavefront) should be
% specified in order to properly apply the map to the wavefront.
% For mirror surface or wavefront error maps, the map values are assumed
% to be in meters, unless the 'nm' or 'microns' switches are used to
% specify the units.  The amplitude map must range from 0 to 1.
% The map will be interpolated to match the current wavefront sampling
% if necessary.
%
% Outputs:
% bm   = wavefront structure with map applied
% map  = 2D output array
%
% Required inputs:
% bm   = the current wavefront structure
% flnm = file name of FITS file, including extension, containing map
% One of the following switches must be specified:
% 'amplitude'         : file contains an amplitude error map
% 'mirror_surface'    : file contains a mirror surface height error map (m)
%                       A positive value indicates a surface point
%                       higher than the mean surface.  The map will be
%                       multiplied by -2 to convert it to a wavefront
%                       map to account for reflection and wavefront
%                       delay (a low region on the surface causes a
%                       positive increase in the phase relative to the
%                       mean).
% 'wavefront'         : file contains a wavefront error map (m) (default)
%
% Optional inputs:
% 'nm'                : map values are in nanometers
% 'magnify'           = spatially magnify the map by this factor from
%                       its default size
% 'microns'           : map values are in microns
% 'multiply'          = multiplies the map amplitude by the specified
%                       factor (default = 1.0)
% 'rotatemap'         = Counter-clockwise rotation of map, after any
%                       resampling and shifting (degrees) (default = 0.0)
% 'sampling'          = sampling of map (m)
% 'xc_map'            = pixel coordinates of map center X
% 'yc_map'            = pixel coordinates of map center Y
% 'xshift'            = amount to shift map in X direction (m) (default = 0.0)
% 'yshift'            = amount to shift map in Y direction (m) (default = 0.0)

% 2005 Feb     jek  created idl routine
% 2014 Aug 18  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed prop_shift_center call to allow odd size arrays
% 2017 Mar 06  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  ampl = 0;             % file does not contain an amplitude map
  cx   = 0.0;           % pixel coordinates of map center X
  cy   = 0.0;           % pixel coordinates of map center Y
  magn = 1.0;           % magnification factor
  mirr = 0;             % file does not contain a mirror map
  mul  = 1.0;           % map amplitude multiplication factor
  nm   = 0;             % map values are not in nanometers
  rot  = 0.0;           % counter-clockwise rotation of map (deg)
  smpl = 0.0;           % sampling of map (m)
  sx   = 0.0;           % amount to shift map in X direction (m)
  sy   = 0.0;           % amount to shift map in Y direction (m)
  um   = 0;             % map values are not in microns
  wave = 1;             % file contains a wavefront error map (default)

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'ampl', 'amplitude'}
        ampl = 1;               % file contains an amplitude map
      case {'cx', 'xc_map'}
        icav = icav + 1;
        cx   = varargin{icav};  % pixel coordinates of map center X
      case {'cy', 'yc_map'}
        icav = icav + 1;
        cy   = varargin{icav};  % pixel coordinates of map center Y
      case {'magn', 'magnify'}
        icav = icav + 1;
        magn = varargin{icav};  % magnification factor
      case {'mul', 'multiply'}
        icav = icav + 1;
        mul  = varargin{icav};  % map amplitude multiplication factor
      case {'mirr', 'mirror_surface'}
        mirr = 1;               % file contains a mirror map
      case {'nm'}
        nm   = 1;               % map values are in nanometers
      case {'rot', 'rotatemap'}
        icav = icav + 1;
        rot  = varargin{icav};  % counter-clockwise rotation of map (deg)
      case {'smpl', 'sampling'}
        icav = icav + 1;
        smpl = varargin{icav};  % sampling of map (m)
      case {'sx', 'xshift'}
        icav = icav + 1;
        sx   = varargin{icav};  % map shift in X direction (m)
      case {'sy', 'yshift'}
        icav = icav + 1;
        sy   = varargin{icav};  % map shift in Y direction (m)
      case {'um', 'microns'}
        um   = 1;               % map values are in microns
      case {'wave', 'wavefront'}
        wave = 1;               % file contains a wavefront error map
      otherwise
        error('prop_errormap: Unknown keyword: %s\n', varargin{icav});
    end
  end

  if exist(flnm) ~= 2           % file does not exist or is not FITS
    error('PROP_ERRORMAP', ...
          'File %s does not exist or is not a FITS file.\n', flnm);
  end

  if ampl == 1 & (nm == 1 | um == 1)
    error('PROP_ERRORMAP', ...
          'Cannot specify units for an amplitude map.\n');
  end

  if cx == 0.0 & cy == 0.0
    info = fitsinfo(flnm);      % FITS file header info
    [lna1, ina1] = ismember('NAXIS1', info.PrimaryData.Keywords(:, 1));
    [lna2, ina2] = ismember('NAXIS2', info.PrimaryData.Keywords(:, 1));
    cx   = fix(info.PrimaryData.Keywords{ina1, 2} / 2) + 1;
    cy   = fix(info.PrimaryData.Keywords{ina2, 2} / 2) + 1;
  end

  if smpl ~= 0.0
    map  = prop_readmap(bm, flnm, ...
      'stx', sx, 'sty', sy, 'mpcx', cx, 'mpcy', cy, 'smpl', smpl);
  else
    map  = prop_readmap(bm, flnm, ...
      'stx', sx, 'sty', sy, 'mpcx', cx, 'mpcy', cy);
  end

  if magn ~= 1.0 | rot ~= 0.0   % apply magnification and/or rotation
    map  = prop_shift_center(map);
    if rot ~= 0.0               % apply rotation
      map  = prop_rotate(map, rot, 'cubic');
    end
    if magn ~= 1.0              % apply magnification
      map  = prop_magnify(map, magn, 'size_out', size(map, 2));
    end
    map  = prop_shift_center(map, 'inv');
  end

  if     nm == 1
    map  = map * 1.0e-9;
  elseif um == 1
    map  = map * 1.0e-6;
  end

  if mul ~= 1.0                 % apply multiplication factor
    map  = map * mul;
  end

  if     ampl == 1              % amplitude error map
    bm.wf = bm.wf .* map;
  elseif mirr == 1              % mirror surface height error map
    bm.wf = bm.wf .* exp(complex(0.0, -4.0 * pi * map / bm.wl));
  elseif wave == 1              % wavefront error map
    bm.wf = bm.wf .* exp(complex(0.0,  2.0 * pi * map / bm.wl));
  else
    error('PROP_ERRORMAP', ...
'Unspecified map type: Use ''mirror_surface'', ''wavefront'', or ''amplitude''');
  end

  map  = prop_shift_center(map);
end                     % function prop_errormap
