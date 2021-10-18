% -------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% function mask = falco_gen_spiral_staircase(inputs)
% 
% Generate a spiral staircase pattern.
%
% FPM must be an array with dimensions MxM, where M = lambdaOverD*(beam diameter)
%
% INPUTS
% ------
% inputs: structure of inputs parameters
%   - inputs.charge: number of 2-pi phase progressions over the 360
%     degrees of the mask). Must be an integer.
%   - inputs.N: width and height of the output array
%   - inputs.phaseScaleFac: Factor to apply uniformly to the phase.
%                           Used to add chromaticity.
%   - inputs.Nsteps: (required for 'staircase') number of steps per
%                    2*pi radians in the staircase
%   - inputs.clocking: (optional) clocking of the phase mask in degrees
%   - inputs.centering: (optional) array centering of the mask
%   - inputs.xOffset: (optional) x-offset of the mask center from the
%                     array center in units of pixels
%   - inputs.yOffset: (optional) y-offset of the mask center from the
%                     array center in units of pixels
%
% OUTPUTS
% -------
% OUT: 2-D, NxN array with complex transmission of the mask 

function mask = falco_gen_spiral_staircase(inputs)

    % Required inputs
    NarrayFinal = inputs.N;
    Nsteps = inputs.Nsteps;
    charge = inputs.charge;
    phaseScaleFac = inputs.phaseScaleFac;

    % OPTIONAL INPUTS
    centering = 'pixel';  %--Default to pixel centering
    xOffset = 0;
    yOffset = 0;
    clocking = 0; % [degrees]
    if(isfield(inputs,'centering')); centering = inputs.centering; end % 'pixel' or 'interpixel'
    if(isfield(inputs, 'xOffset')); xOffset = inputs.xOffset; end % [pixels]
    if(isfield(inputs, 'yOffset')); yOffset = inputs.yOffset; end % [pixels]
    if(isfield(inputs,'clocking')); clocking = inputs.clocking; end % [degrees]

    % Input checks
    Check.scalar_integer(charge);
    
    % minimum radius  of circle centered at (x, y) = (xOffset, yOffset) to contain the full array.
    maxOffset = max(abs(xOffset), abs(yOffset));
    R = 2*NarrayFinal + maxOffset + 1;
    Narray = ceil_even(2*(R+maxOffset) + 2); % Temporary width of array.

    % Set up PROPER
    Dbeam = Narray; %--diameter of beam (normalized to itself)
    dx = Dbeam/Narray;
    Darray = Narray*dx;
    wl_dummy = 1e-6; %--dummy value
    bdf = Dbeam/Darray; %--beam diameter fraction
    wavestruct = prop_begin(Dbeam, wl_dummy, Narray, 'beam_diam_fraction', bdf);

    switch centering
        case {'interpixel'}
            cshift = -dx/2; 
        case {'pixel'}
            cshift = 0;
    end

    deltaThetaStep = 360/charge/Nsteps; % [degrees]
    deltaThetaSector = 360/charge; 
    mask = zeros(Narray);
    for iSector = 1:charge
        for iStep = Nsteps:-1:1
            theta0 = clocking + (iSector-1)*deltaThetaSector - deltaThetaStep;
            xyCorners = [0, R*cosd(theta0), R*cosd(theta0-iStep*deltaThetaStep), 0;
                         0, R*sind(theta0), R*sind(theta0-iStep*deltaThetaStep), 0];

            % Translate
            xyCorners =  xyCorners + cshift;
            xs = xyCorners(1, :) + xOffset;
            ys = xyCorners(2, :) + yOffset;

            % Add next step
            maskNew = (2*pi)/Nsteps*(prop_irregular_polygon(wavestruct, xs, ys));
            mask = mask + maskNew;

        end
    end

    mask = mask - (2*pi)/Nsteps; % Set minimum to zero.
    
    mask = mask * phaseScaleFac;

    mask = pad_crop(mask, NarrayFinal);

end