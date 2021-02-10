%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_wts(bm, dz)
%        bm = prop_wts(bm, dz)
% Propagate a wavefront from a planar reference surface that is
% inside the Rayleigh limit from its focus to a spherical reference
% surface that is inside.  Used by prop_propagate.
% Intended only for use by prop_propagate.  Not a user-callable routine.
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)
% dz   = distance to propagate wavefront (m)

% 2005 Feb     jek  created idl routine
% 2014 May 12  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if abs(dz) <= 1e-12
    return
  end

  global prop_verbose
  if prop_verbose
    fprintf(1, '  WTS: dz:                %10.3f\n', dz);
  end

  if ~strcmp(bm.RefSurf, 'PLANAR')
    fprintf(1, '  WTS: Input reference surface not planar.\n');
    pause
  end
  bm.pz = bm.pz + dz;

  bm   = prop_qphase(bm, dz);

  [ny, nx] = size(bm.wf);
  srn  = sqrt(nx * ny);

  if dz > 0.0          % use forward transform
    bm.wf =  fft2(bm.wf);
% Note that  fft2 includes a normalization of 1
    bm.wf = bm.wf / srn;
  else                  % use inverse transform
    bm.wf = ifft2(bm.wf);
% Note that ifft2 includes a normalization of N^-2
    bm.wf = bm.wf * srn;
  end

  global prop_phase_offset
  if prop_phase_offset
    bm.wf = bm.wf * exp(i * 2 * pi * dz / bm.wl);
  end

  bm.dx = bm.wl * abs(dz) / size(bm.wf, 1) / bm.dx;
  if prop_verbose
    fprintf(1, '  WTS: z:                 %10.3f'  , bm.pz);
    fprintf(1, '       dx:                %10.3e\n', bm.dx);
  end

  bm.RefSurf = 'SPHERI';
end                     % function prop_wts
