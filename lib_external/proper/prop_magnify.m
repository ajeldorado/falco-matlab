%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function aro = prop_magnify(ari, magn, varargin)
%        aro = prop_magnify(ari, magn, varargin)
% function aro = prop_magnify(ari, magn, csrv, nox, quik)
% Resample input array using damped sinc or cubic interpolation
% Calls either an external C routine to do a damped sinc interpolation,
% or uses prop_szoom.m (which uses the same damped sinc interpolation),
% or uses Matlab's interp2 routine to do a cubic interpolation.
%
% Outputs:
% aro  = 2D output array
%
% Required inputs:
% ari  = input array to be magnified
% magn = magnification (e.g., 0.5 = shrink array by factor of 2)
%
% Optional inputs:
% 'amp_conserve'      : If set, the real-valued image is field amplitude
%                       rather than intensity, and to conserve the
%                       resulting intensity, the interpolated image is
%                       divided by the magnification.
% 'conserve'          : If set, the intensity in the image will be
%                       conserved.  If the image is complex, then it is
%                       assumed that the input is an electric field, so
%                       the interpolated result will be divided by the
%                       magnification.  If the image is not complex,
%                       then it is assumed that it is intensity
%                       (modulus-square of the amplitude), so the
%                       interpolated result will be divided by the
%                       square of the magnification.  If the non-complex
%                       image is amplitude and not intensity, specify
%                       'amp-conserve' instead.
% 'quick'             : If set, Matlab's "interp2" function is used
%                       (with "cubic" interpolation) instead of the
%                       more exact, but much slower, sinc interpolator.
%                       Most of the time this will work as well.
% 'size_out'          = dimension of new image (size_out by size_out).
%                       If not specified, the input image dimension is
%                       multiplied by the magnification
% 'quik'  = 1 = Use the external C damped sinc interpolation routine (default)
%         = 2 = Use the prop_szoom damped sinc interpolation routine
%         = 3 = Use the Matlab interp2 cubic interpolation routine

% 2005 Feb     jek  created idl routine
% 2014 Aug 19  gmg  Matlab translation
% 2015 Mar 18  gmg  Added call to interp2
% 2016 Aug 23  gmg  Added call to prop_szoom.m
% 2017 Mar 06  gmg  Revised for keyword/value for optional inputs
% 2018 May 23  jek  Fixed bug that crashed routine if size_out not specified
% 2019 May 28  jek  Fixed incorrect call to fix() 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  ampc = 0;     % image is not field amplitude
  csrv = 0;     % do not conserve intensity
  nox  = 0;     % dimension of new image (pixels)
  if exist('prop_szoom_c') == 3
    quik = 1;   % Use the external C damped sinc interpolation routine
  else
    quik = 2;   % Use the Matlab Proper function prop_szoom.m
  end

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'ampc', 'amp_conserve'}
        ampc = 1;               % conserve intensity
      case {'csrv', 'conserve'}
        csrv = 1;               % conserve intensity
      case {'nox', 'size_out'}
        icav = icav + 1;
        nox  = varargin{icav};  % dimension of new image
      case {'quick'}
        quik = 3;               % use interp2 cubic interpolation
      case {'quik'}
        icav = icav + 1;
        quik = varargin{icav};
      otherwise
        error('prop_magnify: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  if nargin < 2
    error('Proper:PROP_MAGNIFY', ...
          'Must specify input array and magnification.\n');
  end

  if nox == 0
    noy = fix(size(ari,1) .* magn);      		% size of output array
    nox = fix(size(ari,2) .* magn);      		% size of output array
  else
    noy  = nox;
  end

  if(quik == 3)                 % use Matlab interp2 cubic interpolation
                                % works for non-square arrays
    [ny, nx] = size(ari);       % size of input array
    cx   = fix(nx / 2) + 1;     % coordinates of ari center x (pixels)
    cy   = fix(ny / 2) + 1;     % coordinates of ari center y (pixels)
    sx   = fix(nox / 2) + 1;    % coordinates of aro center x (pixels)
    sy   = fix(noy / 2) + 1;    % coordinates of aro center y (pixels)
    [aox, aoy] = meshgrid(([1:nox]-sx) / magn + cx, ([1:noy]-sy) / magn + cy);
    if isreal(ari) == 0         % ari is complex
      aror = zeros(noy, nox);
      aror = interp2(real(ari), aox, aoy, 'cubic');
      aroi = zeros(noy, nox);
      aroi = interp2(imag(ari), aox, aoy, 'cubic');
      aro  = complex(aror, aroi);
    else                        % ari is not complex
      aro  = zeros(noy, nox);
      aro  = interp2(ari, aox, aoy, 'cubic');
    end
  elseif(quik == 2)             % use prop_szoom damped sinc interpolation
                                % works for non-square arrays
    if isreal(ari) == 0         % ari is complex
      aror = zeros(noy, nox);
      aror = prop_szoom(real(ari), magn, 'nox', nox, 'noy', noy);
      aroi = zeros(noy, nox);
      aroi = prop_szoom(imag(ari), magn, 'nox', nox, 'noy', noy);
      aro  = complex(aror, aroi);
    else                        % ari is not complex
      aro  = zeros(noy, nox);
      aro  = prop_szoom(ari, magn, 'nox', nox, 'noy', noy);
    end
  else                          % use external C damped sinc interpolation
                                % assumes square arrays
    if isreal(ari) == 0         % ari is complex
      aror = zeros(nox, nox);
      aror = prop_szoom_c(real(ari), size(ari, 2), nox, magn);
      aroi = zeros(nox, nox);
      aroi = prop_szoom_c(imag(ari), size(ari, 2), nox, magn);
      aro  = complex(aror, aroi);
    else                        % ari is not complex
      aro  = zeros(nox, nox);
      aro  = prop_szoom_c(ari, size(ari, 2), nox, magn);
    end
  end                           % if(quik == 3)

% Conserve intensity?
  if csrv == 1                  % scale array amplitudes
    if isreal(aro) == 1         % array is intensity
      aro  = aro / magn^2;
    else                        % array is electric field
      aro  = aro / magn;
    end
  elseif ampc == 1              % array is amplitude
    aro  = aro / magn;
  end

end                     % function prop_magnify
