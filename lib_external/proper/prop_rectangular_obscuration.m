%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_rectangular_obscuration(bm, sx, sy, varargin)
%        bm = prop_rectangular_obscuration(bm, sx, sy, varargin)
% Multiplies the wavefront in bm by a rectangular obscuration
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)
% sx   = size along x of obscuration (meters unless norm is set)
% sy   = size along y of obscuration (meters unless norm is set)
%
% Optional inputs:
% 'xc'                = center of obscuration relative to the center of
% 'yc'                  the beam (m unless 'norm' is set); default is (0, 0)
% 'rotation'          = counter-clockwise rotation of rectangle about
%                       its center (deg); default is 0.0
% 'norm'              : for sx, sy, xc, yc normalized to beam radius

% 2005 Feb     jek  created idl routine
% 2016 Jan 14  gmg  Matlab translation
% 2017 Feb 28  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  cx   = 0.0;                   % obscuration center to beam center X
  cy   = 0.0;                   % obscuration center to beam center Y
  norm = 0;
  rot  = 0.0;                   % rotation (deg)

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};  % aperture center to beam center X
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};  % aperture center to beam center Y
      case {'norm'}
        norm = 1;
      case {'rot', 'rotation'}
        icav = icav + 1;
        rot  = varargin{icav};  % rotation (degrees)
      otherwise
        error('prop_rectangular_obscuration: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  if norm == 1
  	bm.wf = bm.wf .* prop_shift_center(prop_rectangle(bm, sx, sy, 'cx', cx, 'cy', cy, 'rotation', rot, 'norm', 'dark'), 'inv');
  else
  	bm.wf = bm.wf .* prop_shift_center(prop_rectangle(bm, sx, sy, 'cx', cx, 'cy', cy, 'rotation', rot, 'dark'), 'inv');
  end
end                     % function prop_rectangular_obscuration
