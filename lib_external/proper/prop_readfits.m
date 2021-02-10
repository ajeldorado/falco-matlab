%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [fa] = prop_readfits(flnm)
%        [fa] = prop_readfits(flnm)
% Read in an array from a Flexible Image Transport System (FITS) file
%
% Outputs:
% fa   = output array
%
% Required inputs:
% flnm = file name of FITS file
%
% Intended for internal use by Proper routines.
% Users should call either prop_errormap or prop_psd_errormap.

% 2005 Feb     jek  created idl routine
% 2014 Jul 17  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed fftshift to ifftshift to allow for odd size arrays
% 2016 Jan 15  gmg  Derived from prop_readmap.m
% 2017 Mar 10  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  info = fitsinfo(flnm);        % FITS file header info
  fa   = fitsread(flnm);        % FITS file array

end                     % function prop_readfits
