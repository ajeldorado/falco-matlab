%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function aimo = prop_pixellate(aimi, sami, samo, npo)
%        aimo = prop_pixellate(aimi, sami, samo, npo)
% Integrate a sampled Point Spread Function onto detector pixels.
% This routine takes as input a sampled PSF and integrates it over
% pixels of a specified size.  This is done by convolving the Fourier
% transform of the input PSF with a sinc function representing the
% transfer function of an idealized square pixel and transforming back.
% This result then represents the PSF integrated onto detector-sized
% pixels with the same sampling as the PSF.  The result is interpolated
% to get detector-sized pixels at detector-pixel spacing.
%
% Outputs:
% aimo = 2D array image integrated over square detector pixels
%
% Required inputs:
% aimi = 2D array image containing Point Spread Function
% sami = sampling of aimi (meters / pixel)
% samo = sampling of detector pixels (meters / pixel)
%
% Optional inputs:
% npo  = number of pixels across output image

% 2005 Feb     jek  created idl routine
% 2016 Aug 25  gmg  Matlab translation
% 2017 Mar 21  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  [npiy, npix] = size(aimi);    % number of pixels across input image
  icix = fix(npix / 2) + 1;     % index of center of input image x
  iciy = fix(npiy / 2) + 1;     % index of center of input image y

% Compute pixel transfer function (Modulation Transfer Function)

  magn = sami / samo;           % magnification
  vx   = ifftshift([1 : npix] - icix) / (icix - 1) / 2d0 / magn;
  vy   = ifftshift([1 : npiy] - iciy) / (iciy - 1) / 2d0 / magn;
  [ax, ay] = meshgrid(vx, vy);
  pmtf = sinc(ax) .* sinc(ay);  % pixel modulation transfer function

% Convolve image with detector pixel
  amtf = pmtf .* fft2(ifftshift(aimi));
  cimc = fftshift(abs(ifft2(amtf)) / magn^2);

% Image is integrated over pixels, but has original sampling;
% now resample to pixel sampling
  if (nargin < 4) | (npo < 1)
    npo  = fix(magn * npix);
  end
  aimo = prop_magnify(cimc, magn, 'nox', npo);

end                     % function prop_pixellate
