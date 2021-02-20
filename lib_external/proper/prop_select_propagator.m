%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [bm, dzw] = prop_select_propagator(bm, dz)
%        [bm, dzw] = prop_select_propagator(bm, dz)
% Used by prop_propagate to decide in which propagation regime
% the next surface will be (to decide which propagation method to use
% (spherical-to-planar, planar-to-spherical, or planar-to-planar)).
%
% Outputs:
% bm   = beam structure (output)
% dzw  = distance to new waist position from new position z (m)
%
% Required inputs:
% bm   = beam structure (input)
% dz   = distance to new position z from current position z (m)

% 2005 Feb     jek  created idl routine
% 2014 May 09  gmg  Matlab translation
% 2017 Mar 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  global RayFact

  dzw  = bm.w0_pz - bm.pz;
  newz = bm.pz + dz;

  if abs(bm.w0_pz - newz) < RayFact * bm.zRay
    TypeNew  = 'INSIDE_';
  else
    TypeNew  = 'OUTSIDE';
  end

  bm.PropType = [bm.TypeOld '_to_' TypeNew];
  bm.TypeOld  = TypeNew;
end                     % function prop_select_propagator
