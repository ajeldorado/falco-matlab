% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Object containing the properties of a DM with a square grid of
% surface-perpendicular actuators.

classdef DM
    
    properties
        Nact double {mustBePositive} = 48; % Number of actuators across the width of the square DM grid
        orientation string = 'rot0'; % DM Orienation. Options: rot0, rot90, rot180, rot270, flipxrot0, flipxrot90, flipxrot180, flipxrot270
        xtilt double = 0; % Angle of rotation about the x-axis [degrees] to use in prop_dm.m. For foreshortening the DM surface
        ytilt double = 0; % Angle of rotation about the y-axis [degrees] to use in prop_dm.m. For foreshortening the DM surface
        zrot double = 0; % Clocking of the DM suface [degrees] to use in prop_dm.m.
        xc double = 23.5 % x-center location of DM surface [actuator widths] as given to prop_dm.m. Centered would be at (Nact/2 - 1/2) because of 0-based indexing.
        yc double = 23.5 % x-center location of DM surface [actuator widths] as given to prop_dm.m. Centered would be at (Nact/2 - 1/2) because of 0-based indexing.
        
        dead double = []; % Vector of linear indices of all dead actuators (those stuck at 0V absolute). This should stay fixed.
        pinned double = []; % Vector of linear indices of all actuators to pin at their upper or lower bounds. This is dynamic.
        Vpinned double = []; % Vector of relative voltage commands of pinned/railed actuators, relative to the biasMap.
        tied double = zeros(0, 2); % Indices of tied actuator pairs sharing the same voltage. Two indices per row.
        comovingGroups cell = {}; % Cell array with each index containing a vector of linear indices for actuators that move together. The vectors can be any length.

        Vmin double = 0; % Min allowed absolute voltage command
        Vmax double = 100; % Max allowed absolute voltage command
        dVnbrLat double {mustBeNonnegative} = 100; % max voltage difference allowed between laterally-adjacent DM actuators
        dVnbrDiag double {mustBeNonnegative} = 100; % max voltage difference allowed between diagonally-adjacent DM actuators
        marginNbrRule double = 1e-3; % voltage tolerance used when checking neighbor rule and bound limits. Units of volts.
        
        biasMap double = 50*ones(48); % 2-D, Nact-by-Nact map of absolute voltages. FALCO treats this as the zero point. Needed prior to WFSC to allow + and - voltages. Total voltage is biasMap + V
        V double = zeros(48); % 2-D, Nact-by-Nact map of relative voltages from biasMap
        facesheetFlatmap double = 50*ones(48);% 2-D, Nact-by-Nact map of absolute voltages that produces a flat DM surface. Needed for enforcing the neighbor rule. Locations for dead actuators should be 0.
        
        fitType string = 'linear'; % Type of response for displacement vs voltage. Options are 'linear', 'quadratic', and 'fourier2'.
        VtoH double = 5e-9*ones(48); % 2-D, Nact-by-Nact map of actuator gains in meters/volt.
        
        % Voltage discretization
        HminStepMethod string = 'round'; % Method to use when discretizing voltages. Options are 'round', 'floor', 'ceil', and 'fix'.
        HminStep double {mustBeNonnegative} = 0; % Minimum allowed DM step size in units of meters (not volts).
        
        inf_fn string = 'influence_dm5v2.fits'; % File name of influence function to use. Must have the correct header. Options are 'influence_dm5v2.fits' for one type of Xinetics DM, 'influence_BMC_2kDM_400micron_res10.fits' for BMC 2k DM, and 'influence_BMC_kiloDM_300micron_res10_spline.fits' for BMC kiloDM
        inf_sign string = '+'; % Sign to apply to the influence function for a positive voltage poke. Options are '+' and '-'.
        dm_spacing double {mustBePositive} = 0.9906e-3; % Center-to-center actuator spacing as used in the model. Units of meters.
    
    end % properties
    

    methods

        function S = DM(varargin)
            % Initialize the DM object in 3 possible ways:
            %
            % S = DM
            % S = DM(dmStruct) % where dmStruct is a struct
            % S = DM('property', value, 'property', value, ...)
            
            % Single argument is a struct with fields that will become
            % the object properties.
            if length(varargin) == 1 && isstruct(varargin{1})
                dmStruct = varargin{1}; % for convenience
                fnames = fieldnames(dmStruct);
                for ii = 1:length(fnames)
                    if isprop(S, fnames{ii})
                        S.(fnames{ii}) = dmStruct.(fnames{ii});
                    end
                end
            end % if init with dmStruct struct
            
            % Arguments are explicit pairs of property, value.
            if length(varargin) > 1
                % varargin is property, value pairs
                for ii = 1:2:length(varargin)
                    S.(varargin{ii}) = varargin{ii+1};
                end                 
            end 
    
        end
    
    end % methods
   
end
