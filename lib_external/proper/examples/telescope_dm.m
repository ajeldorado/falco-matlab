%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function wavefront = telescope_dm(wavefront, fl_lens, use_errors, use_dm)
%        wavefront = telescope_dm(wavefront, fl_lens, use_errors, use_dm)
%
% Outputs:
% wavefront  = beam structure (output)
%
% Required inputs:
% wavefront  = beam structure (input)
% fl_lens    = focal length (m)
% use_errors : if set, use prop_psd_errormap
% use_dm     : if set, use a deformable mirror

% 2005 Feb     jek  created idl routine
% 2017 Feb 09  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
% 2017 Oct 11  gmg  Revised to match IDL variable names
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
  wavefront = prop_propagate(wavefront, fl_lens, 'snm', 'DM');

  if use_dm == 1
    nact =    49;               % number of DM actuators along one axis
    nact_across_pupil =    47;  % number of DM actuators across pupil
    dm_xc = fix(nact / 2);      % actuator X index at wavefront center
    dm_yc = fix(nact / 2);      % actuator Y index at wavefront center
    d_beam = 2.0d0 * prop_get_beamradius(wavefront);    % beam diameter
    act_spacing = d_beam / nact_across_pupil;           % actuator spacing
    map_spacing = prop_get_sampling(wavefront);         % map sampling

% Have passed through focus, so pupil has rotated 180 degrees;
% need to rotate error map (also need to shift due to the way
% the rotate function operates to recenter map)

    obj_map = circshift(rot90(obj_map, 2), [1, 1]);

% Interpolate map to match number of DM actuators

    dm_map = prop_magnify(obj_map, map_spacing / act_spacing, ...
                'size_out', nact);

% Need to put on opposite pattern;
% convert wavefront error to surface height

    wavefront = prop_dm(wavefront, -dm_map / 2.0d0, ...
                  dm_xc, dm_yc, act_spacing, 'fit');
  end

  wavefront = prop_propagate(wavefront, fl_lens, 'snm', 'coronagraph lens');

end                     % function telescope_dm
