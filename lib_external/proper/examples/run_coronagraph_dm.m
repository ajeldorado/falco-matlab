%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wavefront,sampling]=run_coronagraph_dm(wavelength,grid_size,optval)
%        [wavefront,sampling]=run_coronagraph_dm(wavelength,grid_size,optval)
%
% Outputs:
% wavefront   = 2D wave front array intensity
% sampling    = 2D wave front array sampling distance (m)
%
% Required inputs:
% wavelength  = wavelength (m)
% grid_size   = number of pixels
% optval      = pass values structure

% 2005 Feb     jek  created idl routine
% 2017 Feb 09  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
% 2017 Oct 11  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  diam        =    0.1d0;               % diameter (m)
  f_lens      =   24.0d0 * diam;        % focal length (m)
  beam_ratio  =    0.3d0;               % beam diameter fraction

  wavefront   = prop_begin(diam, wavelength, grid_size, beam_ratio);
  wavefront   = prop_circular_aperture(wavefront, diam / 2.0d0);
  wavefront   = prop_define_entrance(wavefront);

  wavefront   = telescope_dm(wavefront,f_lens,optval.use_errors,optval.use_dm);
  wavefront   = coronagraph(wavefront, f_lens, optval.occulter_type, diam);

  [wavefront, sampling] = prop_end(wavefront);

end                     % function run_coronagraph_dm
