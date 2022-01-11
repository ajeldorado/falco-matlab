% Copyright 2018-2021 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% Downsample a 2-D array by filtering and then using interp2d().
%
% INPUTS
% ------
% arrayIn: 2-D array to downsample
% downSampFac: Factor by which to downsample. Must be >0 & <=1. 
% centering: array centering of the input and output arrays. 
%   Must be 'pixel' or 'interpixel'.
%
% OUTPUTS
% -------
% arrayIn: 2-D array after downsampling

function arrayOut = falco_filtered_downsample(arrayIn, downSampFac, centering)

    % Input checks
    Check.two_dim_array(arrayIn);
    Check.real_positive_scalar(downSampFac);
    if ~strcmpi(centering, 'pixel') && ~strcmpi(centering, 'interpixel')
        error("centering must be 'pixel' or 'interpixel'.")
    end
    
    if downSampFac > 1
        error('This function is for downsampling only.')
        
    elseif downSampFac == 1
        arrayOut = arrayIn;
        
    else
        
        diamIn = 1;
        diamOut = downSampFac;

        % Coordinates at original resolution
        dx0 = 1/diamIn;
        if strcmpi(centering, 'pixel')
            arrayIn = pad_crop(arrayIn, ceil_odd(max(size(arrayIn)))); % pad to odd-sized square array
        elseif strcmpi(centering, 'interpixel')
            arrayIn = pad_crop(arrayIn, ceil_even(max(size(arrayIn)))); % pad to even-sized square array
        end
        N0 = max(size(arrayIn));
        x0 = (-(N0-1)/2:(N0-1)/2)*dx0;
        [X0, Y0] = meshgrid(x0);
        if strcmpi(centering, 'pixel')
            R0 = sqrt(X0.^2 + Y0.^2);
        elseif strcmpi(centering, 'interpixel')
            R0 = sqrt((X0+dx0/2).^2 + (Y0+dx0/2).^2);
        end
        
        % Coordinates at new resolution
        dx1 = 1/diamOut;
        if strcmpi(centering, 'pixel')
            N1 = ceil_odd(N0*dx0/dx1);
        elseif strcmpi(centering, 'interpixel')
            N1 = ceil_even(N0*dx0/dx1);
        end
        x1 = (-(N1-1)/2:(N1-1)/2)*dx1;
        [X1, Y1] = meshgrid(x1);
        
        % Make window
        Window = 0*R0;
        Window(R0 <= dx1/2) = 1;
        Window = Window/sum(sum(Window));
        % figure(1); imagesc(Window); axis xy equal tight; colormap jet; colorbar;
        
        arrayIn = ifftshift(ifft2(fft2(fftshift(Window)).*fft2(fftshift(arrayIn)))); %--To get good grayscale edges, convolve with the correct window before downsampling.
        arrayIn = circshift(arrayIn, [1 1]); %--Undo a centering shift

        arrayOut = interp2(X0, Y0, arrayIn, X1, Y1, 'cubic', 0); %--Downsample by interpolation

    end
    
end
