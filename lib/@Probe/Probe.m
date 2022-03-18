% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Object of parameters for pairwise probing.

classdef Probe
    
    properties
        Npairs = 3; % Number of pair-wise probe pairs to use.
        whichDM = 1; % Which DM to use for probing. 1 or 2.
        xOffset = 0; % Offset of probe center at DM in x [actuators]. Use to avoid obscurations.
        yOffset = 0; % Offset of probe center at DM in y [actuators]. Use to avoid obscurations.
        gainFudge = 1; % empirical fudge factor to make average probe amplitude match desired value.

        radius = 12; % Half-width of the square probed region in the image plane [lambda/D]. (NOTE: Only used for square probes.)
        axis = 'alternate'; % Which axis to have the phase discontinuity along. Values can be 'x', 'y', or 'xy' / 'alt' / 'alternate'. The 'alternate' option causes the bad axis to switch between x and y for each subsequent probe pair.  (NOTE: Only used for square probes.)
        
        width = 12; % Width of rectangular probe in focal plane [lambda/D]. (NOTE: Only used for rectangular probes. radius is used instead for a square probed region)
        xiOffset = 6; % Horizontal offset from star of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. No offset for square probed region.)
        height = 24; % Height of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. radius is used instead for a square probed region)
        etaOffset = 0; % Vertical offset from star of rectangular probe in focal plane [lambda/D].  (NOTE: Only used for rectangular probes. No offset for square probed region.)

   end
   
end
