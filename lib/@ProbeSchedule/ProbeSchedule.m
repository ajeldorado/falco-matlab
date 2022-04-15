% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Object containing the vectors of probing values scheduled at each WFSC iteration.

classdef ProbeSchedule
    
    properties
        xOffsetVec int8 % Vector of x-offsets (one value per WFSC iteration) of the probe center from the DM grid center [actuators]
        yOffsetVec int8 % Vector of x-offsets (one value per WFSC iteration) of the probe center from the DM grid center [actuators]
        % if defined, InormProbeVec overrides the automated probe intensity
        InormProbeVec double {mustBePositive} %= 1e-6*ones(10, 1) % Vector of the desired normalized intensity of the probes at each WFSC iteration 
   end
   
end
