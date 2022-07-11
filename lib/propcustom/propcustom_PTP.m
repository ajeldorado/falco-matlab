% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Near-field Fresnel propagator. It uses the convolution, transfer function
% form of the Fresnel integral.  
%
% This version requires a matrix with an even number of points in each
% dimension.
%
% The FFT code is from Computational Fourier Optics: A MATLAB Tutorial by
% David Voelz in Chapter 5, page 63-64.
%
% INPUTS
% ------
% Ein = input E-field as 2-D square array
% Larray = width of array [meters]
% lambda = wavelength of light [meters]
% deltaz = propagation distance [meters]
%
% OUTPUTS
% -------
% Eout: output E-field the same size as the input E-field array

function Eout = propcustom_PTP(Ein, Larray, lambda, deltaz)
    
    if abs(deltaz) > lambda
        
        [M, N] = size(Ein);  
        dx = Larray/N;
        Ncritical = floor(lambda*abs(deltaz)/dx^2); % Critical sampling

        if M ~= N
            error('propcustom_PTP: Error: propcustom_PTP needs a square input array.');
        elseif N < Ncritical % Verify that the angular spectrum sampling requirement is not violated
            fprintf('Min allowed (ideal) number of points across = %u, and actual number = %u \n', Ncritical, N);
            error('propcustom_PTP: WARNING, NOT ENOUGH POINTS USED!')
        else

            % Make coordinate vectors and maps in the frequency plane (faster than using meshgrid.m)
            fx = [0:(N/2-1),(-N/2:-1)]/(N*dx); % fftshifted frequency-plane coordinates.  dimensions = 1xN
            FX2 = ones(N, 1)*fx.^2; % squared, shifted frequency-plane coordinates. dimensions = NxN
            H = exp(-1j*pi*lambda*deltaz*(FX2 + FX2.')); % transfer function

            U1 = fft2(fftshift(Ein)); % shift, then fft source field
            Eout = ifftshift(ifft2(H.*U1)); % inverse fft, then center observation field    

        end
        
    else
        Eout = Ein;
    end
    
end