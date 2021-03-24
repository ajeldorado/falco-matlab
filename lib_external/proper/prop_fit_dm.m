%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [dmzc, dmsf] = prop_fit_dm(dmz, infk)
%        [dmzc, dmsf] = prop_fit_dm(dmz, infk)
% Determine deformable mirror actuator piston values
% that generate a desired DM surface, accounting for
% the effect of the actuator influence function.
%
% Outputs:
% dmzc = DM actuator positions that create the desired surface
%        when the influence function is included (2D image)
% dmsf = dmzc array convolved with infk
%
% Required inputs:
% dmz  = DM surface to match (2D array, with each element representing
%        the desired surface amplitude for the corresponding actuator)
% infk = influence function kernel (2D image sampled at actuator spacing)
%
% Intended for use by the prop_dm routine and not for general users

% 2005 Feb     jek  created idl routine
% 2014 Jun 25  gmg  Matlab translation
% 2016 Jun 09  gmg  fixed bugs
% 2017 Feb 15  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  e0   = 1.0e6;                         % initial assumed rms error

  dmzc = dmz;                           % initial actuator positions
  last = dmzc;                          % last good dm
  dmsf = conv2(dmzc, infk, 'same');     % 2D convolution
  diff = dmz - dmsf;
  erms = sqrt(sum(sum(diff.^2)));       % root mean square error

  while erms < e0 & (e0 - erms) / erms > 0.01
    last = dmzc;                        % last good dm
    e0   = erms;                        % previous iteration rms error
    dmzc = dmzc + diff;
    dmsf = conv2(dmzc, infk, 'same');   % 2D convolution
    diff = dmz - dmsf;
    if max(max(abs(diff))) < 1.0e-15
      break
    end
    erms = sqrt(sum(sum(diff.^2)));     % root mean square error
  end

% If fit diverged, use the result from the previous iteration
  if erms > e0
    dmzc = last;
  end

end                             % function prop_fit_dm
