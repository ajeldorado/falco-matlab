% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Object containing the vectors of probing values scheduled at each WFSC iteration.

classdef ProbeSchedule
    
    properties
        xOffsetVec double % Vector of x-offsets (one value per WFSC iteration) of the probe center from the ProbeSchedule grid center [actuators]
        yOffsetVec double % Vector of y-offsets (one value per WFSC iteration) of the probe center from the ProbeSchedule grid center [actuators]
        rotationVec double % Vector of the rotation angle to add to the probes at each WFSC iteration [degrees]
        % if defined, InormProbeVec overrides the automated probe intensity
        InormProbeVec double {mustBePositive} % Vector of the desired normalized intensity of the probes at each WFSC iteration 

    end % properties

    
    methods

        function S = ProbeSchedule(varargin)
            % Initialize the ProbeSchedule object in 3 possible ways:
            %
            % S = ProbeSchedule
            % S = ProbeSchedule(schedStruct) % where schedStruct is a struct
            % S = ProbeSchedule('property', value, 'property', value, ...)
            
            % Single argument is a struct with fields that will become
            % the object properties.
            if length(varargin) == 1 && isstruct(varargin{1})
                schedStruct = varargin{1}; % for convenience
                fnames = fieldnames(schedStruct);
                for ii = 1:length(fnames)
                    if isprop(S, fnames{ii})
                        S.(fnames{ii}) = schedStruct.(fnames{ii});
                    end
                end
            end % if init with schedStruct struct
            
            % Arguments are explicit pairs of property, value.
            if length(varargin) > 1
                % varargin is property, value pairs
                for ii = 1:2:length(varargin)
                    S.(varargin{ii}) = varargin{ii+1};
                end                 
            end 
    
        end % init method
    
    end % methods
   
end
