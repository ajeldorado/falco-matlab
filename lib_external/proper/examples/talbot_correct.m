%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wfa, samp] = talbot_correct(wlm, np, ov)
%        [wfa, samp] = talbot_correct(wlm, np, ov)
% The Talbot Effect example from the Proper Manual
%
% Outputs:
% wfa  = 2D wavefront complex array
% samp = pixel sampling (m)
%
% Required inputs:
% wlm  = wavelength (m)
% np   = number of pixels
% ov   = optional values structure
% ov.diam  = diameter (m)
% ov.dist  = propagation distance (m)
% ov.peri  = period (m)

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tal  = 2.0d0 * ov.peri^2 / wlm;       % Talbot length (m)

% Create beam structure with prop_begin
  bm   = prop_begin(ov.diam, wlm, np);

% Create 1-D grating pattern
  pa   =    0.200d0          ;  % pattern amplitude
  vx   = ((1 : np) - fix (np / 2) - 1) * prop_get_sampling(bm);

% Create 2-D amplitude grating pattern
  [gx, gy] = meshgrid(vx, vx);
  grt  = 0.5d0 * (1.0d0 + pa * cos(2.0d0 * pi * gx / ov.peri));

  bm   = prop_multiply(bm, grt);
  bm   = prop_define_entrance(bm);

  if ov.dist <= (0.25d0 * tal)
    bm   = prop_propagate(bm, ov.dist, 'to_plane');
  else
    bm   = prop_propagate(bm, 0.25d0 * tal, 'to_plane');
    pha  = prop_get_phase(bm);          % in radians
    pha  = pha * wlm / (2.0d0 * pi);    % in meters
    bm   = prop_add_phase(bm, -pha);
    bm   = prop_propagate(bm, ov.dist - 0.25d0 * tal, 'to_plane');
  end

  [wfa, samp] = prop_end(bm, 'noabs');

end                     % function talbot
