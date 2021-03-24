%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


% Test script for run_occulter.m

% 2005 Feb     jek  created idl routine
% 2017 Feb 08  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
% 2017 Nov 13  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Define plot parameters
  global ifig                ;  % index of figure
  ifig   =     0             ;
  n      =   512             ;  % number of pixels
  lambda =     0.550d0       ;  % wavelength (um)

  optval.occulter_type = 'SOLID'     ;  % solid occulter
  solid = prop_run('run_occulter', lambda, n, 'passvalue', optval);

  optval.occulter_type = 'GAUSSIAN';    % gaussian occulter
  gaussian = prop_run('run_occulter', lambda, n, 'passvalue', optval);

  optval.occulter_type = '8TH_ORDER';   % prop_8th_order_mask
  eighth_order = prop_run('run_occulter', lambda, n, 'passvalue', optval);

  fprintf(1, 'Maximum speckle flux / stellar flux :\n');
  fprintf(1, '   SOLID                 = %16.7e\n', max(max(solid)));
  fprintf(1, '   GAUSSIAN              = %16.7e\n', max(max(gaussian)));
  fprintf(1, '   8TH_ORDER             = %16.7e\n', max(max(eighth_order)));
