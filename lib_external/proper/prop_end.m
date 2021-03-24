%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [wfa, samp] = prop_end(bm, varargin)
%        [wfa, samp] = prop_end(bm, varargin)
% Set variables needed to properly conclude a propagation run
%
% Outputs:
% wfa     = wavefront array
% samp    = wavefront array sampling distance (m)
%
% Required inputs:
% bm      = beam structure
%
% Optional inputs:
% 'extract'           = returns the idx by idx central part of wavefront
% 'noabs'             : if set, returns complex wavefront
%                       By default, the intensity (modulus squared) of
%                       the field is returned

% 2005 Feb     jek  created idl routine
% 2014 May 07  gmg  Matlab translation
% 2017 Mar 06  gmg  Revised for keyword/value for optional inputs
% 2017 Apr 21  gmg  Removed sampling distance (arcsecons) return value
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  idx  = 0.0;           % size of central part of wavefront (pixels)
  noab = 0;             % return intensity

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'idx', 'extract'}
        icav = icav + 1;
        idx  = varargin{icav};  % size of central part (pixels)
      case {'noab', 'noabs'}
        noab = 1;
      otherwise
        error('prop_end: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  samp    = prop_get_sampling(bm);

  if noab == 1                          % keep complex E-field values
    wfa   = prop_shift_center(bm.wf);
  else                                  % calculate intensities
    wfa   = prop_shift_center(abs(bm.wf).^2);
  end

  if idx ~= 0                   % select central portion of array
    idxn =  ceil((idx - 1) / 2);
    idxp = floor((idx - 1) / 2);
    i2   = floor(size(bm.wf) / 2) + 1;  % indices of center of array
    wfa  = wfa(i2(1) - idxn : i2(1) + idxp, i2(2) - idxn : i2(2) + idxp);
  end
end                     % function prop_end
