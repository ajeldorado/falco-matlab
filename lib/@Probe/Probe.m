% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Object of parameters for pairwise probing.

classdef Probe
    
    properties
        Npairs (1, 1) double {mustBePositive} = 3; % Number of pair-wise probe pairs to use.
        whichDM (1, 1) int16 {mustBePositive} = 1; % Which Probe to use for probing. 1 or 2.
        xOffset (1, 1) double = 0; % x-offset of the probe center from the Probe grid center [actuators]. Use to avoid obscurations.
        yOffset (1, 1) double = 0; % y-offset of the probe center from the Probe grid center [actuators]. Use to avoid obscurations.
        rotation (1, 1) double = 0; % rotation angle applied to the probe command [degrees]
        gainFudge = 1; % empirical fudge factor to make average probe amplitude match desired value.

        radius (1, 1) double {mustBePositive} = 12; % Half-width of the square probed region in the image plane [lambda/D]. (NOTE: Only used for square probes.)
        axis string = 'alternate'; % Which axis to have the phase discontinuity along. Values can be 'x', 'y', or 'xy' / 'alt' / 'alternate'. The 'alternate' option causes the bad axis to switch between x and y for each subsequent probe pair.  (NOTE: Only used for square probes.)
        
        width double {mustBePositive} = 12; % Width of rectangular probe in focal plane [lambda/D]. (NOTE: Only used for rectangular probes. radius is used instead for a square probed region)
        xiOffset double = 6; % Horizontal offset from star of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. No offset for square probed region.)
        height double {mustBePositive} = 24; % Height of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. radius is used instead for a square probed region)
        etaOffset double = 0; % Vertical offset from star of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. No offset for square probed region.)

    end % properties
   
    methods
       
       function S = Probe(varargin)
            % init Probe object
            %
            % S = Probe
            % S = Probe(struct_probe) % where modvar is a struct
            % S = Probe('property', value, 'property', value, ...)
            
            if length(varargin) == 1 && isstruct(varargin{1})
                % input is a modvar struct to initialize properties
                modvar = varargin{1}; % for convenience
                fnames = fieldnames(modvar);
                for ii = 1:length(fnames)
                    if isprop(S, fnames{ii})
                        S.(fnames{ii}) = modvar.(fnames{ii});
                    end
                end
            end % if init with modvar struct
            
            if length(varargin) > 1
                % varargin is property, value pairs
                for ii = 1:2:length(varargin)
                    S.(varargin{ii}) = varargin{ii+1};
                end                 
            end % explicit property, value
           
       end % class init
       
   end % methods
end
