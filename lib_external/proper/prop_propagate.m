%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_propagate(bm, dz, varargin)
%        bm = prop_propagate(bm, dz, varargin)
% Determine which propagator to use to propagate the current wavefront
% by a specified distance and do it.
%
% Outputs:
% bm       = beam structure (output)
%
% Required inputs:
% bm       = beam structure (input)
% dz       = distance to propagate wavefront (m)
%
% Optional inputs:
% 'surface_name'      = name of surface to which to propagate (string)
% 'to_plane'          : if set, propagation is to a plane

% 2005 Feb     jek  created idl routine
% 2014 May 09  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
% 2017 Nov 17  gmg  Fixed bug that sometimes ignored 'to_plane' switch
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  snm  = '';    % name of surface to which to propagate
  noy  = 0;     % number of points in output array Y (pixels)
  tpln = 0;     % 'to_plane' is not set

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'snm', 'surface_name'}
        icav = icav + 1;
        snm  = varargin{icav};  % dimension of new image
      case {'to_plane'}
        tpln = 1;
      otherwise
        error('prop_propagate: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  propcommon

  if print_it
    if strcmp('snm', '')
      fprintf(1, 'Propagating\n');
    else
      fprintf(1, 'Propagating to %s\n', snm);
    end
  end

  [bm, dzw] = prop_select_propagator(bm, dz);
  z2   = bm.pz + dz;

  if tpln
    bm.PropType = [bm.PropType(1 : 11) 'INSIDE_'];
  end

  if prop_verbose
    fprintf(1, '  PROP_PROPAGATOR: propagator_type = %18s\n', bm.PropType);
  end

  if     strcmp(bm.PropType, 'INSIDE__to_INSIDE_')
    bm   = prop_ptp(bm, dz                 );

  elseif strcmp(bm.PropType, 'INSIDE__to_OUTSIDE')
    bm   = prop_ptp(bm, bm.w0_pz - bm.pz   );
    bm   = prop_wts(bm, z2       - bm.w0_pz);

  elseif strcmp(bm.PropType, 'OUTSIDE_to_INSIDE_')
    bm   = prop_stw(bm, bm.w0_pz - bm.pz   );
    bm   = prop_ptp(bm, z2       - bm.w0_pz);

  elseif strcmp(bm.PropType, 'OUTSIDE_to_OUTSIDE')
    bm   = prop_stw(bm, bm.w0_pz - bm.pz   );
    bm   = prop_wts(bm, z2       - bm.w0_pz);
  end

  if print_total_intensity
    tint = sum(sum(abs(bm.wf).^2));
    if strcmp('snm', '')
      fprintf(1, 'Total intensity:          %10.3e\n', tint);
    else
      fprintf(1, 'Total intensity at surface %s: %10.3e\n', snm, tint);
    end
  end

% Save stuff for layout plots

  if do_table
    ActionNum = ActionNum + 1;  % list index
    bmdl(ActionNum) = 2 * prop_get_beamradius(bm);
%   bmdl = list of beam diameters at each lens (m)
    dzl(ActionNum)  = dz;       % list of propagation distances (m)
    saml(ActionNum) = bm.dx;    % list of sampling at each surface (m)
    if strcmp('snm', '')
      snml{ActionNum} = '(Surface)';
    else
      snml{ActionNum} = snm;    % list of surface names
    end
  end
end                     % function prop_propagate
