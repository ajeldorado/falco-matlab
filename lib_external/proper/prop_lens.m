%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_lens(bm, fl, snm)
%        bm = prop_lens(bm, fl, snm)
% Alter the current wavefront as a perfect lens would.
% This routine computes the phase change caused by a perfect lens
% that has a focal length specified by the user.
% A positive focal length corresponds to a convex lens or concave mirror;
% a negative focal length corresponds to a concave lens or convex mirror.
% This routine updates the new beam waist position.
%
% Outputs:
% bm   = wavefront structure (output)
%
% Required inputs:
% bm   = wavefront structure (input)
% fl   = focal length of lens (m)
%
% Optional inputs:
% snm  = surface name (string) (optional)

% 2005 Feb     jek  created idl routine
% 2014 May 08  gmg  Matlab translation
% 2016 Sep 12  gmg  Changed prop_radius2 to prop_radius
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon
  if print_it
    if nargin < 3
      fprintf(1, 'Applying lens\n');
    else
      fprintf(1, 'Applying lens at %s\n', snm);
    end
  end

  bm.zRay = pi * bm.w0^2 / bm.wl;
  wsrf = bm.w0 * sqrt(1.0 + ((bm.pz - bm.w0_pz) / bm.zRay)^2);
% wsrf = beam radius at lens surface (m)

  dzw0 = bm.pz - bm.w0_pz;
% fprintf(1, '  Lens: dzw0:       %16.9e'  , dzw0);
% fprintf(1, '        wsrf:       %16.9e\n', wsrf);

  if dzw0 ~= 0.0        % lens is not at focus or entrance pupil
    gRbo = dzw0 + bm.zRay^2 / dzw0;
    gRboinf  = 0.0;
    if gRbo ~= fl
      gRbm = 1.0 / (1.0 / gRbo - 1.0 / fl);
      gRbminf  = 0.0;   % output beam is spherical
      if prop_verbose
        fprintf(1, '  LENS: Gaussian Rbeam old:%13.6e'  , gRbo);
        fprintf(1, '        R beam:     %20.13e\n', gRbm);
      end
    else
      gRbminf  = 1.0;   % output beam is planar
      if prop_verbose
        fprintf(1, '  LENS: Gaussian Rbeam old:%13.6e'  , gRbo);
        fprintf(1, '        R beam:                 infinite\n');
      end
    end
  else                  % at focus or entrance pupil, input beam is planar
    gRboinf  = 1.0;
    gRbm = -fl;
    gRbminf  = 0.0;     % output beam is spherical
    if prop_verbose
      fprintf(1, '  LENS: Gaussian Rbeam old:     infinite');
      fprintf(1, '        R beam:     %20.13e\n', gRbm);
    end
  end

  if strcmp(bm.TypeOld, 'INSIDE_') | strcmp(bm.RefSurf, 'PLANAR')
    Rbo  = 0.0;
  else
    Rbo  = dzw0;
  end

  if not(gRbminf)
    bm.w0_pz = -gRbm / (1.0 + (bm.wl * gRbm / (pi * wsrf^2))^2) + bm.pz;
    bm.w0 = wsrf / sqrt(1.0 + (pi * wsrf^2 / (bm.wl * gRbm))^2);
  else                  % output beam is planar
    bm.w0_pz = bm.pz;
    bm.w0 = wsrf;
  end
% fprintf(1, '  Lens: w0:         %16.9e'  , bm.w0);
% fprintf(1, '        w0_pz:      %16.9e\n', bm.w0_pz);

% Determine new Rayleigh distance from focus;
% if currently inside this, then output beam is planar

  bm.zRay = pi * bm.w0^2 / bm.wl;
  dzw0 = bm.pz - bm.w0_pz;
% fprintf(1, '  Lens: dzw0:       %16.9e'  , dzw0);
% fprintf(1, '        zRay:       %16.9e\n', bm.zRay);

  if abs(dzw0) < RayFact * bm.zRay
    TypeNew  = 'INSIDE_';
    Rbm  = 0.0;
  else
    TypeNew  = 'OUTSIDE';
    Rbm  = dzw0;
  end

  bm.PropType = [bm.TypeOld '_to_' TypeNew];

% Apply phase changes as needed, but don't apply if the phase
% is going to be similarly altered during propagation.

  if prop_verbose
    fprintf(1, '  LENS: propagator type:        %18s\n', bm.PropType);
    if strcmp(bm.TypeOld, 'INSIDE_')
      sRbo = 'Infinite';
    else
      sRbo = sprintf('%8.3f', Rbo);
    end
    if strcmp(TypeNew, 'INSIDE_')
      sRbm = 'Infinite';
    else
      sRbm = sprintf('%8.3f', Rbm);
    end
    fprintf(1, '  LENS: Rbmold: %8s  ', sRbo);
    fprintf(1, '        R_beam: %8s  ', sRbm);
    fprintf(1, '  focal length: %12.7f\n', fl);
    fprintf(1, '  LENS: beam diam at lens:    %20.13e\n', 2.0 * wsrf);
  end

  if     strcmp(bm.PropType, 'INSIDE__to_INSIDE_')
    phas = 1.0 / fl;
 %  fprintf(1, 'INSIDE__to_INSIDE_  phas:   %16.9e\n', phas);

  elseif strcmp(bm.PropType, 'INSIDE__to_OUTSIDE')
    phas = 1.0 / fl + 1.0 / Rbm;
 %  fprintf(1, 'INSIDE__to_OUTSIDE  phas:   %16.9e  Rbm:    %16.9e\n', phas, Rbm);

  elseif strcmp(bm.PropType, 'OUTSIDE_to_INSIDE_')
    phas = 1.0 / fl - 1.0 / Rbo;
 %  fprintf(1, 'OUTSIDE_to_INSIDE_  phas:   %16.9e  Rbo:    %16.9e\n', phas, Rbo);

  elseif strcmp(bm.PropType, 'OUTSIDE_to_OUTSIDE')
    if     Rbo == 0.0
      phas = 1.0 / fl + 1.0 / Rbm;
    elseif Rbm == 0.0
      phas = 1.0 / fl - 1.0 / Rbo;
    else
      phas = 1.0 / fl - 1.0 / Rbo + 1.0 / Rbm;
      if prop_verbose
        fprintf(1, '  LENS: 1/lens_fl:            %20.13e\n', 1.0 / fl );
        fprintf(1, '  LENS: 1/R_beam_old:         %20.13e\n', 1.0 / Rbo);
        fprintf(1, '  LENS: 1/R_beam:             %20.13e\n', 1.0 / Rbm);
        fprintf(1, '  LENS: lens_phase:           %20.13e\n', phas);
      end
    end
  end

  bm   = prop_add_phase(bm, -(prop_radius(bm).^2) * phas / 2.0);

  if strcmp(TypeNew, 'INSIDE_')
    bm.RefSurf  = 'PLANAR';
  else
    bm.RefSurf  = 'SPHERI';
  end
  bm.TypeOld  = TypeNew;
  bm.fr       = abs(dzw0) / wsrf / 2.0;

% Save stuff for layout plots

  if do_table
    ActionNum = ActionNum + 1;  % list index
    bmdl(ActionNum) = 2 * wsrf; % list of beam diameters at each lens (m)
    efrl(ActionNum) = bm.fr;    % effective focal ratios after each lens
    fll(ActionNum)  = fl;       % list of lens focal lengths (m)
    saml(ActionNum) = bm.dx;    % list of sampling at each surface (m)
    if nargin > 2
      snml{ActionNum} = snm;    % list of surface names
    else
      snml{ActionNum} = '(LENS)';
    end
  end

  if prop_verbose
    fprintf(1, '  LENS: Rayleigh distance:    %20.13e\n', bm.zRay);
  end
end                     % function prop_lens
