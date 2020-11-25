function arrayOut = offcenter_crop(arrayIn, centerRow, centerCol, nRowOut, nColOut)
% """
% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% Crop a 2-D array to be centered at the specified pixel.
% 
% This function crops a 2-D array about the center pixel specified by
% centerRow and centerCol. The input array can be
% rectangular with even or odd side lengths. The output will is rectangular
% with dimensions nRowOut, nColOut. If the output array includes
% regions outside the input array, those pixels are included
% and set to zero. If the specified cropping region is fully outside the
% input array, then the output is all zeros.
% 
% The center pixel of an odd-sized array is the array
% center, and the center pixel of an even-sized array follows the FFT
% center pixel convention.
% 
% Parameters
% ----------
% arrayIn : array_like
%     2-D input array
% nRowOut, nColOut : int
%     Height and width of the 2-D output array in pixels. Must be integers.
% centerRow, centerCol : float or int
%     Indices of the pixel to be used as the output array's center.
%     Floating point values are rounded to the nearest integer.
%     Convention in this function is that y is the first axis (rows).
%     Values can be negative and/or lie outside the input array.
% 
% Returns
% -------
% arrayOut : numpy ndarray
%     2-D square array
% 
% Notes
% -----
% All alignment units are in pixels.
% """

[nRowIn, nColIn] = size(arrayIn);
centerRow = round(centerRow);
centerCol = round(centerCol);

% Compute how much to pad the array in y (if any)
rowPadPre = 0;
rowPadPost = 0;
if ceil(-nRowOut/2) + centerRow < 1
    rowPadPre = abs(ceil(-nRowOut/2) + centerRow);
end
if ceil(nRowOut/2.) + centerRow > nRowIn
    rowPadPost = ceil(nRowOut/2.) + centerRow - nRowIn;
end
rowPad = 1 + max([rowPadPre, rowPadPost]);

% Compute how much to pad the array in x (if any)
colPadPre = 0;
colPadPost = 0;
if ceil(-nColOut/2.) + centerCol < 1
    colPadPre = abs(ceil(-nColOut/2) + centerCol);
end
if ceil(nColOut/2.) + centerCol > nColIn
    colPadPost = ceil(nColOut/2) + centerCol - nColIn;
end
colPad = 1 + max([colPadPre, colPadPost]);

arrayPadded = pad_crop(arrayIn, [nRowIn+2*rowPad, nColIn+2*colPad]);

centerRow = centerRow + rowPad - 1;
centerCol = centerCol + colPad - 1;

% Buffer needed to keep output array correct size
if mod(nRowOut, 2) == 1
    rowBuffer = 1;
elseif mod(nRowOut, 2) == 0
    rowBuffer = 0;
end
if mod(nColOut, 2) == 1
    colBuffer = 1;
elseif mod(nColOut, 2) == 0
    colBuffer = 0;
end

arrayOut = arrayPadded(...
   centerRow-floor(nRowOut/2)+1:centerRow+floor(nRowOut/2)+rowBuffer,...
   centerCol-floor(nColOut/2)+1:centerCol+floor(nColOut/2)+colBuffer);

end
