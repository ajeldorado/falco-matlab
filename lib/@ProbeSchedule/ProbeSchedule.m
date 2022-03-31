% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Object containing the vectors of probing values scheduled at each WFSC iteration.

classdef ProbeSchedule
    
    properties
        xOffsetVec double % Vector of x-offsets (one value per WFSC iteration) of the probe center from the DM grid center [actuators]
        yOffsetVec double % Vector of y-offsets (one value per WFSC iteration) of the probe center from the DM grid center [actuators]
        rotationVec double % Vector of the rotation angle to add to the probes at each WFSC iteration [degrees]
        InormProbeVec double {mustBePositive} % Vector of the desired normalized intensity of the probes at each WFSC iteration 
   end
   
end
