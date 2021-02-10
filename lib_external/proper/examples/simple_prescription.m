%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wfai, samp] = simple_prescription(wlm, np)
%        [wfai, samp] = simple_prescription(wlm, np)
%
% Outputs:
% wfai = 2D wave front array intensity
% samp = 2D wave front array sampling distance (m)
%
% Required inputs:
% wlm  = wavelength (m)
% np   = number of pixels

% 2005 Feb     jek  created idl routine
% 2017 Apr 21  gmg  Matlab translation
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Default parameters
  bdf  =    0.500d0          ;  % beam diameter fraction
  diam =    1.000d0          ;  % telescope diameter (m)
  fr   =   15.000d0          ;  % focal ratio
  flm  = diam * fr           ;  % focal length (m)

  bm   = prop_begin(diam, wlm, np, bdf);

% Create circular entrance aperture with 0.5 m radius
  bm   = prop_circular_aperture(bm, diam / 2.0d0);
  bm   = prop_define_entrance(bm);

  bm   = prop_lens(bm, flm);
  bm   = prop_propagate(bm, flm);

% Calculate 2D wave front array intensity and sampling distance (m)
  [wfai, samp] = prop_end(bm);

end                     % function simple_prescription
