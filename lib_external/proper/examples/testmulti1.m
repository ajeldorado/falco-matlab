%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


% Demo script for testmulti1.m

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
% 2017 Nov 14  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  final_sampling =    1.500d-6; % final sampling (m)
  fn1            = 'output_testmulti1.fits';
  gridsize       = 1024       ; % number of pixels
  npsf           =  256       ; % number of pixels, Point Spread Functions
  nlambda        =    9       ; % number of wavelengths
  lambda_min     =    0.500d0 ; % wavelength minimum (um)
  lambda_max     =    0.700d0 ; % wavelength maximum (um)

% Generate array of wavelengths (um)

  lambda_um=lambda_min+(lambda_max-lambda_min)*[0:(nlambda-1)]/(nlambda-1);

% Create Deformable Mirror pattern (a couple of 0.1 micron pokes)

  dm         = zeros(48, 48);
  dm(21, 21) =    0.200d-6  ;
  dm(16, 26) =    0.200d-6  ;
  use_dm     =    1         ;   % multi_example flag: use deformable mirror
  optval     = struct('use_dm', use_dm, 'dm', dm);

% Generate monochromatic fields in parallel

  [fields,sampling]=prop_run_multi('multi_example',lambda_um,gridsize,'passvalue',optval);

% Resample fields to same scale, convert to Point Spread Functions

  psfs = zeros(npsf, npsf, nlambda);
  for i = 1 : nlambda
    mag   = sampling(i) / final_sampling;
    field = prop_magnify(fields(:, :, i), mag, 'size_out', npsf, 'conserve');
    psfs(1 : npsf, 1 : npsf, i) = abs(field).^2;
  end

% Add PSFs together

  psf  = sum(psfs, 3) / nlambda;
  fprintf(1, '     total: %14.7e      mean: %14.7e       max: %14.7e\n', ...
          sum(sum(psf)), mean(mean(psf)), max(max(psf)));
  prop_writemap(psf, fn1);
