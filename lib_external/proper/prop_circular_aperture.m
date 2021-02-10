%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_circular_aperture(bm, ra, varargin)
%        bm = prop_circular_aperture(bm, ra, varargin)
% Multiplies the wavefront in bm by a circular aperture
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)
% ra   = radius of aperture (meters unless 'norm' is set)
%
% Optional inputs:
% 'xc'                = center coordinate X (meters unless 'norm')
% 'yc'                = center coordinate Y (meters unless 'norm')
% 'norm'              : for ra, cx, cy normalized to beam radius

% 2005 Feb     jek  created idl routine
% 2014 Jun 09  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed prop_shift_center call to allow odd size arrays
% 2017 Mar 06  gmg  Revised for keyword/value for optional inputs
% 2020 Jan 30  jek  Fixed NORM in call to prop_ellipse
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  cx   = 0.0;                   % aperture center to beam center X
  cy   = 0.0;                   % aperture center to beam center Y
  norm = 0;

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
      otherwise
        error('prop_circular_aperture: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  if ( norm == 1 )
  	bm.wf = bm.wf .* prop_shift_center(prop_ellipse(bm, ra, ra, 'cx', cx, 'cy', cy, 'norm'), 'inv');
  else
  	bm.wf = bm.wf .* prop_shift_center(prop_ellipse(bm, ra, ra, 'cx', cx, 'cy', cy), 'inv');
  end

end                     % function prop_circular_aperture
