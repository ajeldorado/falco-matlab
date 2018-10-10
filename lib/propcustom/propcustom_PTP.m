% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = prop_PTP(Ein,Larray,lambda,dz) 
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
% Larray = width of array (meters)
% lambda = wavelength of light (meters)
% dz = propagation distance (meters)
%
% OUTPUTS:
%  Eout: output E-field as square array
%
% Modified on 2018-07-27 by A.J. Riggs not to use meshgrid.
% Modified on 2018-04-26 by A.J. Riggs to correct the sampling requirement.
% Written by A.J. on April 30, 2013. 

function Eout = propcustom_PTP(Ein,Larray,lambda,dz)
    % propagation - transfer function approach
    % assumes same x and y side lengths
    % uniform sampling
    % Ein - source plane field
    % Larray - source and observation plane side length
    % lambda - wavelength
    % dz - propagation distance
    % u2 - observation plane field
    [M,N] = size(Ein);  
    dx = Larray/N;
    Ncritical = floor(lambda*abs(dz)/dx^2); %--Critical sampling
    
    if(M~=N) %--Check if input array is square
        error('propcustom_PTP: Error: propcustom_PTP needs a square input array.');
    elseif( N < Ncritical) %--Verify that the angular spectrum sampling requirement is not violated
        fprintf('Min allowed (ideal) number of points across = %u, and actual number = %u \n', Ncritical, N);
        error('propcustom_PTP: WARNING, NOT ENOUGH POINTS USED!')
    else
%         fx = ( -M/2:1:(M/2-1) )/Larray;   % dx = Larray/M; fx=-1/(2*dx):1/Larray:1/(2*dx)-1/Larray;    % freq coords
%         [FX,FY]=meshgrid(fx,fx);
%         H = fftshift( exp(-1i*pi*lambda*dz*(FX.^2+FY.^2)) ); % transfer function

        %--Make coordinate vectors and maps in the frequency plane (faster than using meshgrid.m)
        fx = [0:(N/2-1),(-N/2:-1)]/(N*dx);    % fftshifted frequency-plane coordinates.  dimensions = 1xN
        FX2 = ones(N,1)*fx.^2; % squared, shifted frequency-plane coordinates. dimensions = NxN
        H = exp(-1i*pi*lambda*dz*(FX2+FX2.')); % transfer function

        U1=fft2(fftshift(Ein));      % shift, then fft source field
        Eout=ifftshift(ifft2(H.*U1));    % inv fft, then center observation field    

    end
end
