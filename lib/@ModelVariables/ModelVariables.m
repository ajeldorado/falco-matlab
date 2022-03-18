% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Object of model variables passed to model_compact() or model_full() when
% computing the electric field.

classdef ModelVariables
    
    properties
        
        sbpIndex = 1 % Index of which subband to use
        wpsbpIndex = 1 % Index of which wavelength to use within the subband
        starIndex = 1 % Index of which star to use. Typically is 1 since there is usually just 1 star.
        zernIndex = 1 % Index of which Noll Zernike mode to include as amplitude at the entrance pupil. Typically just want 1 (for piston).
        whichSource = 'star' % Which source. 'star' to put at star location, or 'offaxis' for a temporary offset from the star center.
        x_offset = 0  % x-offset of the source from the nominal star center in lambda0/D. Only used when whichSource = 'offaxis'.
        y_offset = 0  % y-offset of the source from the nominal star center in lambda0/D. Only used when whichSource = 'offaxis'.
        lambda {mustBePositive} % Wavelength at which to compute the E-field. Units of meters. If specified, overrides the values given by sbpIndex and wpsbpIndex. Used only in standalone calls to the model.
        ebpIndex {mustBePositive} % Used only for ZWFS
   end

end
