function Ipixelated = fourierPixelate(I,Npix,varargin)
%Ipixellated = fourierPixelate(I,Npix,varargin)
%   Given intensity distribution I, simulates the measurement by a 2D 
%   detector array with (Npix x Npix) pixels
%   Only works for square arrays with even numbers of rows and columns 
%
%   Inputs:
%       I - Input array (square, even dimensions)
%       Npix - Number of simulated pixels across I 
%       varargin(1) - centering offset [rows,cols] in units of pixels 
%       
%   Outputs:
%       Ipixelated - pixelated image (Npix x Npix)

    [Nrows,Ncols] = size(I);

    if( mod(Nrows,2)~=0 || mod(Ncols,2)~=0 || Nrows~=Ncols  )
        error('fourierPixellate: Input array must be even and square');
    end

    if(numel(varargin)>0)
        centerPixOffset = varargin{1};
    else
        centerPixOffset = [0 0];
    end
    
    Narr = Nrows; 
    sampsPerPix = Narr/Npix; 
    
    %%-- Low pass filter 
    [X,Y] = meshgrid(-Narr/2:Narr/2-1); 
    a = Narr/sampsPerPix;
    ftCameraPixel = fftshift(sinc(X/a).*sinc(Y/a));
    I_lpf = ifft2(fft2(I).*ftCameraPixel);
    
    %%-- Resample at pixel centers 
    [Xpix, Ypix] = meshgrid(-Npix/2:Npix/2-1);

    Ipixelated = interp2(X, Y, I_lpf,...
        (Xpix-centerPixOffset(2))*sampsPerPix,...
        (Ypix-centerPixOffset(1))*sampsPerPix);
    
end

