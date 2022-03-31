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
        whichDM (1, 1) int16 {mustBePositive} = 1; % Which DM to use for probing. 1 or 2.
        xOffset (1, 1) int16 = 0; % x-offset of the probe center from the DM grid center [actuators]. Use to avoid obscurations.
        yOffset (1, 1) int16 = 0; % y-offset of the probe center from the DM grid center [actuators]. Use to avoid obscurations.
        gainFudge = 1; % empirical fudge factor to make average probe amplitude match desired value.

        radius (1, 1) double {mustBePositive} = 12; % Half-width of the square probed region in the image plane [lambda/D]. (NOTE: Only used for square probes.)
        axis string = 'alternate'; % Which axis to have the phase discontinuity along. Values can be 'x', 'y', or 'xy' / 'alt' / 'alternate'. The 'alternate' option causes the bad axis to switch between x and y for each subsequent probe pair.  (NOTE: Only used for square probes.)
        
        width double {mustBePositive} = 12; % Width of rectangular probe in focal plane [lambda/D]. (NOTE: Only used for rectangular probes. radius is used instead for a square probed region)
        xiOffset double = 6; % Horizontal offset from star of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. No offset for square probed region.)
        height double {mustBePositive} = 24; % Height of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. radius is used instead for a square probed region)
        etaOffset double = 0; % Vertical offset from star of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. No offset for square probed region.)

   end
   
end
