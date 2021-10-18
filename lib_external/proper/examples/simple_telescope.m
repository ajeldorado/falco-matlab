%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wfa, samp] = simple_telescope(wlm, np)
%        [wfa, samp] = simple_telescope(wlm, np)
% A Simple Telescope example from the Proper Manual
%
% Outputs:
% wfa  = 2D wavefront intensity array
% samp = pixel sampling (m)
%
% Required inputs:
% wlm  = wavelength (m)
% np   = number of pixels

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dobj =    0.060d0          ;  % objective diameter (m)
  fl1  =   15.000d0 * dobj   ;  % focal length objective (m)
  fl2  =    0.021d0          ;  % focal length eyepiece (m)
  fl3  =    0.022d0          ;  % focal length human eye (m)
  br   =    0.500d0          ;  % initial beam width/grid width

% Create beam structure with prop_begin
  bm   = prop_begin(dobj, wlm, np, br);

  bm   = prop_circular_aperture(bm, dobj / 2.0d0);
  bm   = prop_define_entrance(bm);

  bm   = prop_lens(bm, fl1, 'objective');

  bm   = prop_propagate(bm, fl1 + fl2, 'surface_name', 'eyepiece');
  bm   = prop_lens(bm, fl2, 'eyepiece');

% Calculate exit pupil distance
  epd  = fl2 / (1.0d0 - fl2 / (fl1 + fl2));
  bm   = prop_propagate(bm, epd, 'surface_name', 'exit pupil at eye lens');
  bm   = prop_lens(bm, fl3, 'eye');

  bm   = prop_propagate(bm, fl3, 'surface_name', 'retina');

% Extract 2D array of intensity and pixel sampling from beam structure
  [wfa, samp] = prop_end(bm);

end                     % function simple_telescope
