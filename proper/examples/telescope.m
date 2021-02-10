%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function wavefront = telescope(wavefront, fl_lens, use_errors)
%        wavefront = telescope(wavefront, fl_lens, use_errors)
%
% Outputs:
% wavefront  = beam structure (output)
%
% Required inputs:
% wavefront  = beam structure (input)
% fl_lens    = focal length (m)
% use_errors : if set, use prop_psd_errormap

% 2005 Feb     jek  created idl routine
% 2017 Feb 08  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
% 2017 Nov 20  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if use_errors == 1
    rms_error  =    10.0d-09  ; % RMS wavefront error
    c_freq     =    15.0d0    ; % correlation frequency (cycles / m)
    high_power =     3.0d0    ; % high frequency falloff
    flnm = 'telescope_obj.fits';
    [wavefront, obj_map] = prop_psd_errormap(wavefront, rms_error, ...
                             c_freq, high_power, 'file', flnm, 'rms');
  end

  wavefront = prop_lens(wavefront, fl_lens, 'objective');

% Propagate through focus to the pupil

  wavefront = prop_propagate(wavefront, fl_lens * 2.0d0, ...
                'snm', 'telescope pupil imaging lens');
  wavefront = prop_lens(wavefront, fl_lens, 'telescope pupil imaging lens');

% Propagate to a deformable mirror (no actual DM here)

  wavefront = prop_propagate(wavefront, fl_lens, 'snm', 'DM');
  wavefront = prop_propagate(wavefront, fl_lens, 'snm', 'coronagraph lens');

end                     % function telescope
