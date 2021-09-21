% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% function OUT = propcustom_mft_PtoFtoP(IN, FPM, charge, apRad, inVal, outVal, useGPU ) 
%
% Generate a 
%
% FPM must be an array with dimensions MxM, where M = lambdaOverD*(beam diameter)
%
% INPUTS
% ------
% inputs: structure of inputs parameters
%   - inputs.type: type of mask. Valid options are 'vortex', 'cos',
%   'sectors', and 'staircase.'
%   - inputs.charge: number of 2-pi phase progressions over the 360
%     degrees of the mask)
%   - inputs.N: width and height of the output array
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

function mask = falco_gen_azimuthal_phase_mask(inputs)

    % Required inputs
    maskType = inputs.type; 
    charge = inputs.charge; 
    N = inputs.N; 

    % OPTIONAL INPUTS
    centering = 'pixel';  %--Default to pixel centering
    xOffset = 0;
    yOffset = 0;
    clocking = 0; % [degrees]
    if(isfield(inputs,'centering')); centering = inputs.centering; end % 'pixel' or 'interpixel'
    if(isfield(inputs, 'xOffset')); xOffset = inputs.xOffset; end % [pixels]
    if(isfield(inputs, 'yOffset')); yOffset = inputs.yOffset; end % [pixels]
    if(isfield(inputs,'clocking')); clocking = inputs.clocking; end % [degrees]

    % make coordinate system 
    if strcmpi(centering, 'pixel')
        xs = -N/2:(N/2-1);
    elseif strcmpi(centering, 'interpixel')
        xs = -(N-1)/2:(N-1)/2;
    else
        error('Invalid option for inputs.centering')
    end
    [X, Y] = meshgrid(xs);
    
    X = X - xOffset;
    Y = Y - yOffset;
    
    % rotate all points before computing THETA
    xyAll = [reshape(X, [1, N*N]); reshape(Y, [1, N*N])];
    rotMat = [cosd(clocking), sind(clocking); -sind(clocking), cosd(clocking)];
    xyAll = rotMat * xyAll;
    X = reshape(xyAll(1, :), [N, N]);
    Y = reshape(xyAll(2, :), [N, N]);

    [THETA, ~] = cart2pol(X, Y);

    % make mask 
    switch lower(maskType)
        case 'vortex'
            mask = exp(1j*charge*THETA);

        case{'cos'}
            z_m = besselzero(0,1);
            z_m = z_m(end);
            mask = exp(1j*z_m*cos(charge*THETA));

        case 'sectors'
            mask = exp(1j*pi/2*sign(cos(charge*THETA)));

        case 'staircase'
            Nsteps = inputs.Nsteps;
            mask = exp(1j*ceil(mod((THETA+pi)/(2*pi)*charge, 1)*Nsteps)/Nsteps*2*pi);

        otherwise
            error('%s not a valid option for inputs.type.', inputs.type)
    end

end

