%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_elliptical_aperture(bm, rx, ry, varargin)
%        bm = prop_elliptical_aperture(bm, rx, ry, varargin)
% Multiplies the wavefront in bm by an elliptical clear aperture
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)
% rx   = X radius of aperture (meters unless 'norm' is set)
% ry   = Y radius of aperture (meters unless 'norm' is set)
%
% Optional inputs:
% 'xc'                = center coordinate X (meters unless 'norm')
% 'yc'                = center coordinate Y (meters unless 'norm')
% 'norm'              : for rx, ry, xc, yc normalized to beam radius
% 'rotation'          = angle to rotate in degrees
% 2005 Feb     jek  created idl routine
% 2016 Apr 14  gmg  Matlab translation
% 2017 Mar 01  gmg  Revised for keyword/value for optional inputs
% 2020 Jan 30  jek  Fixed 'norm' in call to prop_ellipse, added rotation, fixed yc/cy parameter
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  cx   = 0.0;                   % aperture center to beam center X
  cy   = 0.0;                   % aperture center to beam center Y
  norm = 0;
  rot = 0;			% rotation angle in degrees
  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};  % aperture center to beam center X
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};  % aperture center to beam center X
      case {'rot', 'rotation'}
        icav = icav + 1;
        rot   = varargin{icav};  % aperture center to beam center Y
      case {'norm'}
        norm = 1;
      otherwise
        error('prop_elliptical_aperture: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  if ( norm == 1 )
	bm.wf = bm.wf .* prop_shift_center(prop_ellipse(bm, rx, ry, 'cx', cx, 'cy', cy, 'rotation', rot, 'norm'), 'inv');
  else
	bm.wf = bm.wf .* prop_shift_center(prop_ellipse(bm, rx, ry, 'cx', cx, 'cy', cy, 'rotation', rot ), 'inv');
  end

end                     % function prop_elliptical_aperture
