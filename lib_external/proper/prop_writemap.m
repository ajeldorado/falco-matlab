%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function prop_writemap(map, flnm, varargin)
%        prop_writemap(map, flnm, varargin)
% Write a wavefront, surface, or amplitude aberration map to a
% Flexible Image Transport System (FITS) file.
%
% Required inputs:
% map  = 2D input array
% flnm = file name of FITS file, including extension
%
% Optional inputs:
% 'amplitude'         : set if amplitude map
% 'mirror'            : set if mirror surface (not wavefront) map in meters
% 'wavefront'         : set if wavefront map in meters (default)
% 'radius_pix'        = Specifies the beam radius in the map in pixels.
%                       If specified, the value of sampling (if provided) is
%                       ignored.  When this file is read by prop_errormap,
%                       the map will be resampled as necessary to match
%                       the sampling of the beam.
% 'sampling'          = Map sampling (m); ignored if radius_pix is specified.

% 2005 Feb     jek  created idl routine
% 2014 Aug 05  gmg  Matlab translation
% 2017 Feb 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if exist(flnm) == 2           % file already exists
    error('Proper:PROP_WRITEMAP', ...
          'File %s already exists.\n', flnm);
  end

  mptp = 'wavefront' ;          % wavefront map
  pixr =     0.0d0   ;          % beam radius (pixels)
  smpl =     0.0d0   ;          % spacing (m)

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'ampl', 'amplitude'}
        mptp = 'amplitude';     % amplitude map
      case {'mirr', 'mirror'}
        mptp = 'mirror'   ;     % mirror map
      case {'wave', 'wavefront'}
        mptp = 'wavefront';     % wavefront map
      case {'pixr', 'radius_pix'}
        icav = icav + 1;
        pixr = varargin{icav};  % beam radius (pixels)
      case {'smpl', 'sampling'}
        icav = icav + 1;
        smpl = varargin{icav};  % map sampling (m)
      otherwise
        error('prop_writemap: Unknown keyword: %s\n', varargin{icav});
    end
  end

  import matlab.io.*;
  fptr = fits.createFile(flnm);

  [ny, nx] = size(map);         % number of rows and columns in map
  fits.createImg(fptr, 'double_img', [ny nx]);
  fits.deleteKey(fptr, 'EXTEND');
  fits.writeDate(fptr);
  fits.writeImg(fptr, map);

  fits.writeKey(fptr, 'MAPTYPE', mptp, ' error map type');

  fits.writeKey(fptr, 'X_UNIT', 'meters', ' X & Y units');
  if ~strcmp(mptp, 'amplitude')
    fits.writeKey(fptr, 'Z_UNIT', 'meters', ' Error units');
  end

  if pixr > 0
    fits.writeKey(fptr, 'RADPIX', pixr, ' beam radius in pixels');
  elseif smpl > 0
    fits.writeKey(fptr, 'PIXSIZE', smpl, ' spacing in meters', 11);
  end

  icpx = int16(floor(nx / 2) + 1);      % index of center pixel x
  icpy = int16(floor(ny / 2) + 1);      % index of center pixel y
  fits.writeKey(fptr, 'XC_PIX', icpx, ' Center X pixel coordinate');
  fits.writeKey(fptr, 'YC_PIX', icpy, ' Center Y pixel coordinate');

  fits.closeFile(fptr);
end                     % function prop_writemap
