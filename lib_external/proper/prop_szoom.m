%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function arro = prop_szoom(arri, magn, varargin)
%        arro = prop_szoom(arri, magn, varargin)
% Perform image magnification about the center using damped sinc
% (Lanczos) interpolation.
% This is intended to be called from prop_magnify.
%
% Outputs:
% arro = 2D output array
%
% Required inputs:
% arri = array to be magnified (maybe real or complex)
% magn = magnification (>1 magnifies image, <1 shrinks)
%
% Optional inputs:
% 'nox'               = number of points in output array X
% 'noy'               = number of points in output array Y (default: nox)

% 2006 Sep     jek  created idl routine
% 2014 Mar     jek  modified to support complex values
% 2016 Aug 15  gmg  Matlab translation; works correctly for non-square arrays
% 2017 Mar 10  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  nox  = 0;     % number of points in output array X (pixels)
  noy  = 0;     % number of points in output array Y (pixels)

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'nox'}
        icav = icav + 1;
        nox  = varargin{icav};  % dimension of new image
      case {'noy'}
        icav = icav + 1;
        noy  = varargin{icav};  % dimension of new image
      otherwise
        error('prop_szoom: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  if nargin < 2
    error('Proper:PROP_SZOOM', ...
          'Must specify input array and magnification.\n');
  end

  fsp  = 6;                     % filter size parameter (input pixels)
  ks   = 1 + 2 * fsp;           % kernel size (input pixels)
  [niy, nix] = size(arri);      % input array size
  cplx = any(imag(arri(:)));    % complex input array
  if nox <= 0
    nox  = fix(nix * magn) - ks;% output array size x
  end
  if noy <= 0
    noy  = nox;
  end
  iicx = fix(nix / 2) + 1;      % index of  input array center x (pixels)
  iocx = fix(nox / 2) + 1;      % index of output array center x (pixels)
  iicy = fix(niy / 2) + 1;      % index of  input array center y (pixels)
  iocy = fix(noy / 2) + 1;      % index of output array center y (pixels)

% Precompute damped sinc kernel table
  vk   = [0 : ks - 1] - fix(ks / 2);    % kernel vector
  vp   = ([0 : nox - 1] - fix(nox / 2)) / magn;
  vp   = vp - round(vp);                % phase vector
  [ak, ap] = meshgrid(vk, vp);
  at   = ak - ap;
  mask = abs(at) <= fsp;
  at   = at * pi;
  au   = ones(nox, ks);
  au(at ~= 0) =  sin(at(at ~= 0)      ) ./ (at(at ~= 0)      ) ...
              .* sin(at(at ~= 0) / fsp) ./ (at(at ~= 0) / fsp);
  tblx = au .* mask;

  if noy ~= nox
    vk   = [0 : ks - 1] - fix(ks / 2);  % kernel vector
    vp   = ([0 : noy - 1] - fix(noy / 2)) / magn;
    vp   = vp - round(vp);              % phase vector
    [ak, ap] = meshgrid(vk, vp);
    at   = ak - ap;
    mask = abs(at) <= fsp;
    at   = at * pi;
    au   = ones(noy, ks);
    au(at ~= 0) =  sin(at(at ~= 0)      ) ./ (at(at ~= 0)      ) ...
                .* sin(at(at ~= 0) / fsp) ./ (at(at ~= 0) / fsp);
    tbly = au .* mask;
  else
    tbly = tblx;
  end

  if cplx
    arro = complex(zeros(noy, nox), zeros(noy, nox));
  else
    arro = zeros(noy, nox);
  end

% Apply kernels
  for ioy = 1 : noy                     % index of output array y
    iiy  = round((ioy-iocy)/magn)+iicy; % index of  input array y
    iiy1 = iiy - fsp;                   % index of  input array y1
    iiy2 = iiy + fsp;                   % index of  input array y2
    if (iiy1 < 1) | (iiy2 > niy)
      continue
    end
    [atx, aty] = meshgrid([1 : nix], tbly(ioy, :));
    strp = arri(iiy1 : iiy2, :) .* aty;
% Note: in Matlab, code is the same for complex and real cases!
    for iox = 1 : nox                           % index of output array x
      iix  = round((iox-iocx)/magn)+iicx;       % index of  strp  array x
      iix1 = iix - fsp;                         % index of  strp  array x1
      iix2 = iix + fsp;                         % index of  strp  array x2
      if(iix1 < 1 | iix2 > nix)
        continue
      end
      ur   = sum(strp(:, iix1 : iix2), 1) .* tblx(iox, :);
      arro(ioy, iox) = sum(ur);
    end                                 % for iox = 1 : nox
  end                           % for ioy = 1 : noy
end                     % function prop_szoom
