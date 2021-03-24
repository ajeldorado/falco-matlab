%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [map] = prop_readmap(bm, flnm, varargin)
%        [map] = prop_readmap(bm, flnm, varargin)
% Read in a surface, wavefront, or amplitude error map from a
% Flexible Image Transport System (FITS) file, scaling if necessary.
%
% Outputs:
% map  = output array
%
% Required inputs:
% bm   = beam structure
% flnm = file name of FITS file
%
% Optional inputs:
% 'xshift'            = amount to shift map in x (m)
% 'yshift'            = amount to shift map in y (m)
% 'xc_map'            = map center pixel x (default is fix(mnx / 2) + 1)
% 'yc_map'            = map center pixel y (default is fix(mny / 2) + 1)
% 'sampling'          = sampling of map (m) (will override any sampling
%                       specified in the file header; must be specified
%                       if header does not specify sampling using the
%                       PIXSIZE value).
%        Note: if the header value RADPIX is specified (the radius of
%        the beam in the map in pixels), then this will override any
%        other sampling specifiers, either in the header or using smpl.
%
% Intended for internal use by Proper routines.
% Users should call either prop_errormap or prop_psd_errormap.

% 2005 Feb     jek  created idl routine
% 2014 Jul 17  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed fftshift to ifftshift to allow for odd size arrays
% 2017 Feb 23  gmg  Revised for keyword/value for optional inputs
% 2017 Aug 22  jek  Fixed error in call to prop_resamplemap found by Arthur Newman and reported by Gerry Rafanelli (Raytheon)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  stx  = 0.0;                   % Set default map shifts to 0.0
  sty  = 0.0;

  info = fitsinfo(flnm);        % FITS file header info
  map  = fitsread(flnm);        % FITS file array

  [mny, mnx] = size(map);       % Set default map center indices
  mpcx = fix(mnx / 2.0) + 1;
  mpcy = fix(mny / 2.0) + 1;

  smpl = 0.0;

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'stx', 'xshift'}
        icav = icav + 1;
        stx  = varargin{icav};  % map shift X (m)
      case {'sty', 'yshift'}
        icav = icav + 1;
        sty  = varargin{icav};  % map shift Y (m)
      case {'mpcx', 'xc_map'}
        icav = icav + 1;
        mpcx = varargin{icav};  % center pixel X (pixels)
      case {'mpcy', 'yc_map'}
        icav = icav + 1;
        mpcy = varargin{icav};  % center pixel Y (pixels)
      case {'smpl', 'sampling'}
        icav = icav + 1;
        smpl = varargin{icav};  % map sampling (m)
      otherwise
        error('prop_readmap: Unknown keyword: %s\n', varargin{icav});
    end
  end

% If the radius of the beam (in pixels) in the map is specified in
% the header (RADPIX value), this will override any other sampling
% specifications, either from the header (PIXSIZE value) or procedure
% call (smpl value).  If smpl is set, it overrides PIXSIZE.

  [lrdp, irdp] = ismember('RADPIX' , info.PrimaryData.Keywords(:, 1));
  [lpxs, ipxs] = ismember('PIXSIZE', info.PrimaryData.Keywords(:, 1));
  if lrdp > 0                   % radpix is defined
    pxsz = prop_get_beamradius(bm) / info.PrimaryData.Keywords{irdp, 2};
  elseif smpl > 0               % smpl is defined
    pxsz = smpl;
  elseif lpxs > 0               % pixsize is defined
    pxsz = info.PrimaryData.Keywords{ipxs, 2};
  else
    error('Proper:PROP_READMAP', ...
          'No pixel scale specified in %s header', flnm);
  end

% Resample map to current wavefront grid spacing
  map  = prop_resamplemap(bm, map, pxsz, mpcx, mpcy, stx, sty);

% Shift center of map to (1, 1)
  map  = ifftshift(map);
end                     % function prop_readmap
