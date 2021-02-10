%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_begin(diam, wl, np, varargin)
%        bm = prop_begin(diam, wl, np, varargin)
% Initialize variables for PROPER routines.  This routine must be called
% before any other PROPER routines in order to initialize required
% variables.
% Creates a structure with a complex array of electric field magnitudes
% = wfa (if present) or = (1.0, 0.0i) otherwise
%
% Outputs:
% bm      = beam structure created by this routine
%
% Required inputs:
% diam    = diameter of beam (m)
% wl      = wavelength (m)
% np      = number of points in each dimension of np x np array
%
% Optional inputs:
% 'beam_diam_fraction'= beam diameter fraction (default = 0.5)

% 2005 Feb     jek  created idl routine 
% 2014 Apr 18  gmg  Matlab translation
% 2017 Feb 14  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% beam structure:
% bm.diam     = initial beam diameter (m)
% bm.dx       = grid sampling (m / pixel)
% bm.fr       = focal ratio
% bm.PropType = 'INSIDE__to_OUTSIDE', 'OUTSIDE_to_INSIDE_', 'INSIDE__to_INSIDE_'
% bm.pz       = beam position z (m)
% bm.Rbeam    = beam radius of curvature (m)
% bm.RbeamInf = beam starts out with infinite curvature radius
% bm.RefSurf  = reference surface type: 'PLANAR' or 'SPHERI'
% bm.TypeOld  = beam location 'INSIDE_', 'OUTSIDE' (inside, outside beam waist)
% bm.w0       = beam waist radius (m)
% bm.w0_pz    = beam waist position z (m)
% bm.wf       = initialized wavefront array
% bm.wl       = wavelength (m)
% bm.zRay     = Rayleigh distance for current beam (m)

  propcommon
  npa         = np;     % number of points across array
  OldOPD      = 0.0;
  RayFact     = 2.0;    % wavefront switches from planar to spherical
                        % at RayFact * bm.zRay
  bdf         = 0.5;    % beam diameter fraction

  if nargin == 4
    bdf  = varargin{1};
  elseif nargin > 4
    icav = 0;                   % index in cell array varargin
    while icav < size(varargin, 2)
      icav = icav + 1;
      switch lower(varargin{icav})
        case {'bdf', 'beam_diam_fraction'}
          icav = icav + 1;
          bdf  = varargin{icav};
        otherwise
          error('prop_begin: Unknown keyword: %s\n', varargin{icav});
      end
    end
  end

  bm.diam     = diam;
  bm.dx       = diam / bdf / np;
  bm.fr       = 1.0d9;          % default is collimated beam
  bm.PropType = 'INSIDE__to_INSIDE_';
  bm.pz       = 0.0;
  bm.Rbeam    = 0.0;
  bm.RbeamInf = 1;              % beam starts with infinite curvature radius
  bm.RefSurf  = 'PLANAR';
  bm.TypeOld  = 'INSIDE_';
  bm.w0       = diam / 2.0;
  bm.w0_pz    = 0.0;
  bm.wf       = complex(ones(np, np), zeros(np, np));
  bm.wl       = wl;
  bm.zRay     = pi * bm.w0^2 / wl;

  antialias_subsampling = int16(11);	% default subsampling in each dimension of shape edges 

  nlist = 1500;
  if do_table
    ActionNum = 0;              % list index
    bmdl = zeros(nlist, 1);     % list of beam diameters at each lens (m)
    dzl  = zeros(nlist, 1);     % list of propagation distances (m)
    efrl =  ones(nlist, 1);     % effective focal ratios after each lens
    fll  = zeros(nlist, 1);     % list of lens focal lengths (m)
    saml = zeros(nlist, 1);     % list of sampling at each surface (m)
    snml = cellstr(char(zeros(nlist, 1)));      % list of surface names
  end

end                     % function prop_begin
