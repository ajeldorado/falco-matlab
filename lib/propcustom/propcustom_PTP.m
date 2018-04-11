% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = prop_PTP(Ein,Darray,lambda,dz) 
%
% This is a Fresnel propagation function. It uses the convolution, transfer function
% form of the Fresnel integral.  
%
% This version requires a matrix with an even number of points in each
% dimension.
%
% The FFT code is from Computational Fourier Optics: A MATLAB Tutorial by
% David Voelz in Chapter 5, page 63-64.
%
%
% REQUIRED INPUTS: 
% Ein = input E-field as 2-D square array
% Darray = width of array (meters)
% lambda = wavelength of light (meters)
% dz = propagation distance (meters)
%
% OUTPUTS:
%  Eout: output E-field as square array
%
% Written by A.J. on April 30, 2013. 

function Eout = propcustom_PTP(Ein,Darray,lambda,dz)
    % propagation - transfer function approach
    % assumes same x and y side lengths
    % uniform sampling
    % Ein - source plane field
    % Darray - source and observation plane side length
    % lambda - wavelength
    % dz - propagation distance
    % u2 - observation plane field
    [M,N] = size(Ein);  
    if(M~=N) %--Check if input array is square
        disp('Error: propcustom_PTP needs a square input array.');
        return;
    elseif(round(Darray^2/lambda/abs(dz))<M) %--Verify that the angular spectrum sampling requirement is not violated
        disp(' WARNING: TOO MANY POINTS USED in function propcustom_PTP')
        fprintf('Max allowed (ideal) number of points across = %u, and actual number = %u \n', floor(Darray^2/lambda/abs(dz)) , M);
    else
        fx = ( -M/2:1:(M/2-1) )/Darray;   % dx = Darray/M; fx=-1/(2*dx):1/Darray:1/(2*dx)-1/Darray;    % freq coords
        [FX,FY]=meshgrid(fx,fx);

        H = fftshift( exp(-1i*pi*lambda*dz*(FX.^2+FY.^2)) ); % transfer function
        U1=fft2(fftshift(Ein));      % shift, then fft source field
        Eout=ifftshift(ifft2(H.*U1));    % inv fft, then center observation field    

    end
end
