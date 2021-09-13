% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% surfaceToFitSmoothed = smooth_beyond_edges(surfaceToFit, windowWidth, nIterSmooth)
%
% Smooth a surface map beyond its edges to make it easier to fit a DM commands to.
%
% INPUTS
% ------
% surfaceToFit : structure of all model parameters
% mask : boolean mask of which pixels correspond to the original surface
%     map
% windowWidth : Width of the convolution kernel in pixels. Rounded up to
%     next odd number.
% nIterSmooth : number of iterations to 
%
% OUTPUTS
% -------
% surfaceToFitSmoothed : 2-D map with edges extended via smoothing

function surfaceToFitSmoothed = smooth_beyond_edges(surfaceToFit, mask, windowWidth, nIterSmooth)

    windowWidth = ceil_odd(windowWidth); % pixels

    kernel = ones(windowWidth);
    kernel = kernel/sum(kernel(:));

    surfaceToFitSmoothed = surfaceToFit; % initialize
    for ii = 1:nIterSmooth
        surfaceToFitSmoothed = conv2(surfaceToFitSmoothed, kernel, 'same');
        surfaceToFitSmoothed(mask) = surfaceToFit(mask); % replace with original values
    end
end
