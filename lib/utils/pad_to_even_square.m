% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Zero-pad a 2-D array to make it square-shaped with even-length sides.
%
% INPUTS
% ------
% arrayIn: 2-D input array
%
% OUTPUTS
% -------
%  arrayOut: 2-D square output array with even side lengths. 

function arrayOut = pad_to_even_square(arrayIn)

    if length(size(arrayIn)) ~= 2
        error('arrayIn must be a 2-D array')
    end
        
    arrayOut = pad_crop(arrayIn, ceil_even(max(size(arrayIn))));
   
end
