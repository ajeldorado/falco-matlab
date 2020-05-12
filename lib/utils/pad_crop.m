% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to pad or crop the input matrix to the desired size.
%
% If both sizes of a dimension are even or odd, the array will be centered
% in that dimension.  If the input size is even and the output size is odd,
% the smaller array will be shifted one element toward the start of the
% array.  If the input is odd and the output is even, the smaller array will
% be shifted one element toward the end of the array.  This sequence is
% intended to ensure transitivity, so several pad_crop calls can be
% chained in no particular order without changing the final result.
%
% Originally written in Python by Brian Kern and Eric Cady.
%
%--INPUTS:
% arrayIn: rectangular matrix to be padded or cropped. Can have even or odd
% side lengths. Can be up to a 4-D array, but only the first two dimensions
% are padded.
% nDes: Number of points desired across the array. If an int, the output is an nDes x nDes array.
%       If a vector, output is a nDes(2) x nDes(1) array. 
%
%--OPTIONAL INPUTS:
% 'EXTRAPVAL', VALUE : The keyword EXTRAPVAL indicates that the array, if padded anywhere, will be
%      padded with the value VALUE instead of the default of 0.
%
%--OUTPUTS:
% arrayOut: rectangular padded or cropped version of the input matrix


function arrayOut = pad_crop(arrayIn, nDes, varargin)

    % Set default values of input parameters
    extrapval = 0; %--Value to use for extrapolated points

    %--Enable different extrapolation values by using varargin
    icav = 0; % index in cell array varargin
    while icav < size(varargin, 2)
        icav = icav + 1;
        switch lower(varargin{icav})
          case {'extrapval'}
            icav = icav + 1;
            extrapval   = varargin{icav};  %--Value to use for extrapolated points
          otherwise
            error('pad_crop: Unknown keyword: %s\n', ...
              varargin{icav});
        end
    end

    nInX = size(arrayIn,2);
    nInY = size(arrayIn,1);
    
    if(length(nDes) == 1)
        nOutX = nDes;
        nOutY = nDes;
    elseif(length(nDes) == 2)
        nOutX = nDes(2);
        nOutY = nDes(1);
    end

    arrayOut = extrapval * ones(nOutY, nOutX, size(arrayIn,3), size(arrayIn,4));

    xneg = min([floor(nInX/2), floor(nOutX/2)]);
    xpos = min([nInX - floor(nInX/2), nOutX - floor(nOutX/2)]);
    yneg = min([floor(nInY/2), floor(nOutY/2)]);
    ypos = min([nInY - floor(nInY/2), nOutY - floor(nOutY/2)]);

    ySliceIn = floor(nInY/2)-yneg+1 : floor(nInY/2)+ypos;
    xSliceIn = floor(nInX/2)-xneg+1 : floor(nInX/2)+xpos;
    
    ySliceOut = floor(nOutY/2)-yneg+1 : floor(nOutY/2)+ypos;
    xSliceOut = floor(nOutX/2)-xneg+1 : floor(nOutX/2)+xpos;

    arrayOut(ySliceOut, xSliceOut, :, :) = arrayIn(ySliceIn, xSliceIn, :, :);
    
end %--END OF FUNCTION
