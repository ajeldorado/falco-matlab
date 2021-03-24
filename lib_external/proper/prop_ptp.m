%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_ptp(bm, dz)
%        bm = prop_ptp(bm, dz)
% This routine is used by prop_propagate to propagate a planar input
% wavefront over some distance to produce a planar output wavefront.
% This occurs when both the start and end point are both within the
% Rayleigh distance of focus.
% Intended only for use by prop_propagate.  Not a user-callable routine.
%
% Outputs:
% bm       = beam structure (output)
%
% Required inputs:
% bm       = beam structure (input)
% dz       = distance to propagate wavefront (m)

% 2005 Feb     jek  created idl routine
% 2013 Oct     jek  speed up phase term calculation
% 2014 May 12  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed fftshift to ifftshift to allow for odd size arrays
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if abs(dz) <= 1e-12
    return
  end

  global prop_verbose
  if prop_verbose
    fprintf(1, '  PTP: dz:                %10.3f\n', dz);
  end

  if ~strcmp(bm.RefSurf, 'PLANAR')
    fprintf(1, '  PTP: Input reference surface not planar.\n');
    pause
  end
  bm.pz = bm.pz + dz;

  [ny, nx] = size(bm.wf);
  srn  = sqrt(nx * ny);
  bm.wf =  fft2(bm.wf);
% Note that  fft2 includes a normalization of 1
  bm.wf = bm.wf / srn;

  ix2 = [-floor(nx / 2) : floor((nx - 1) / 2)].^2;
  iy2 = [-floor(ny / 2) : floor((ny - 1) / 2)].^2;
  [ax2, ay2]  = meshgrid(ix2, iy2);
  ri2 = ifftshift((ax2 + ay2) / nx / ny / prop_get_sampling(bm)^2);
  bm.wf = bm.wf .* exp(-i * pi * bm.wl * dz * ri2);

  bm.wf = ifft2(bm.wf);
% Note that ifft2 includes a normalization of N^-2
  bm.wf = bm.wf * srn;

  global prop_phase_offset
  if prop_phase_offset
    bm.wf = bm.wf * exp(i * 2.0 * pi * dz / bm.wl);
  end

  if prop_verbose
    fprintf(1, '  PTP: z:                 %10.3f'  , bm.pz);
    fprintf(1, '       dx:                %10.3e\n', bm.dx);
  end

end                     % function prop_ptp
