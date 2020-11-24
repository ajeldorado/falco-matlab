function recentered_image = offcenter_crop(arrayIn, pixel_count_across, output_center_x, output_center_y)
% """
% Crop a 2-D array to be centered at the specified pixel.
% 
% This function crops a 2-D array about the center pixel specified by
% output_center_x and output_center_y. The input array can be
% rectangular with even or odd side lengths. The output will be a
% square of side length pixel_count_across. If the output array has
% regions outside the original 2-D array, those pixels are included
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
% pixel_count_across : int
%     Width of the 2-D output array in pixels
% output_center_x, output_center_y : float or int
%     Indices of the pixel to be used as the output array's center.
%     Floating point values are rounded to the nearest integer.
%     Convention in this function is that y is the first axis.
%     Values can be negative and/or lie outside the input array.
% 
% Returns
% -------
% recentered_image : numpy ndarray
%     2-D square array
% 
% Notes
% -----
% All alignment units are in pixels.
% """

[pixel_count_y, pixel_count_x] = size(arrayIn);
x_center = round(output_center_x);
y_center = round(output_center_y);

% Compute how much to pad the array in x (if any)
x_pad_pre = 0;
x_pad_post = 0;
if ceil(-pixel_count_across/2.) + x_center < 1
    x_pad_pre = abs(ceil(-pixel_count_across/2) + x_center);
end
if ceil(pixel_count_across/2.) + x_center > pixel_count_x
    x_pad_post = ceil(pixel_count_across/2) + x_center - pixel_count_x;
end
x_pad = 1 + max([x_pad_pre, x_pad_post]);

% Compute how much to pad the array in y (if any)
y_pad_pre = 0;
y_pad_post = 0;
if ceil(-pixel_count_across/2) + y_center < 1
    y_pad_pre = abs(ceil(-pixel_count_across/2) + y_center);
end
if ceil(pixel_count_across/2.) + y_center > pixel_count_y
    y_pad_post = ceil(pixel_count_across/2.) + y_center - pixel_count_y;
end
y_pad = 1 + max([y_pad_pre, y_pad_post]);

padded_image = pad_crop(arrayIn, [pixel_count_y+2*y_pad, pixel_count_x+2*x_pad]);

x_center = x_center + x_pad - 1;
y_center = y_center + y_pad - 1;

% Buffer needed to keep output array correct size
if mod(pixel_count_across, 2) == 1
    buffer = 1;
elseif mod(pixel_count_across, 2) == 0
    buffer = 0;
end

halfAcross = floor(pixel_count_across/2);
recentered_image = padded_image(...
   y_center-halfAcross+1:y_center+halfAcross+buffer,...
   x_center-halfAcross+1:x_center+halfAcross+buffer);

end
