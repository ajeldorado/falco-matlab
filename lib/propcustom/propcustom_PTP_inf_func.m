% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = propcustom_PTP_inf_func(Ein,Larray,lambda,dz,actPitch,propMethod)
%
% This is a Fresnel propagation function. It uses the convolution, transfer function
% form of the Fresnel integral.  
%
% This version requires a matrix with an even number of points in each
% dimension.
%
% The FFT code is from Computational Fourier Optics: A MATLAB Tutorial by
% David Voelz in Chapter 5, page 63-64. It has been modified not to use
% meshgrid.m for speed.
%
%
% REQUIRED INPUTS: 
% Ein = input E-field as 2-D square array
% Larray = width of array [meters]
% lambda = wavelength of light [meters]
% dz = propagation distance [meters]
% actPitch = inter-actuator pitch of the deformable mirror [meters]
% propMethod = whether to use MFTs or FFTs in the propagation ('mft' or
% 'fft')
%
% OUTPUTS:
%  Eout: output E-field as square array
%
% Modified on 2019-03-11 by A.J. Riggs to be for influence function
% propagation only.
% Modified on 2018-07-27 by A.J. Riggs not to use meshgrid.
% Modified on 2018-04-26 by A.J. Riggs to correct the sampling requirement.
% Written by A.J. on April 30, 2013. 

function Eout = propcustom_PTP_inf_func(Ein,Larray,lambda,dz,actPitch,propMethod)
    
    [M,N] = size(Ein);  
    dx = Larray/N;
    Ncritical = floor(lambda*abs(dz)/dx^2); %--Critical sampling
    
    %--Error checking of inputs
    if(M~=N) %--Check if input array is square
        error('propcustom_PTP: Error: propcustom_PTP needs a square input array.');
    elseif( N < Ncritical) %--Verify that the angular spectrum sampling requirement is not violated
        fprintf('Min allowed (ideal) number of points across = %u, and actual number = %u \n', Ncritical, N);
        error('propcustom_PTP: WARNING, NOT ENOUGH POINTS USED!')
    end
    
    switch lower(propMethod)
        case{'fft'}
            
            %--Make coordinate vectors and maps in the frequency plane (faster than using meshgrid.m)
            fx = [0:(N/2-1),(-N/2:-1)]/(N*dx); % fftshifted frequency-plane coordinates.  dimensions = 1xN
            FX2 = ones(N,1)*fx.^2; % squared, shifted frequency-plane coordinates. dimensions = NxN
            H = exp(-1i*pi*lambda*dz*(FX2+FX2.')); % transfer function

            U1=fft2(fftshift(Ein)); % shift, then fft source field
            Eout=ifftshift(ifft2(H.*U1)); % inv fft, then center observation field    
        
        case{'mft'}
            
            NactPerMat = Larray/actPitch;
            Nxi = ceil_even(NactPerMat*4);
            centering = 'pixel';
            f = 1; %--Dummy, normalized focal length [meters]
            dxi = (f*lambda/Larray);
            
            %--Make coordinate vectors and maps in the frequency plane (faster than using meshgrid.m)
            fx = (-Nxi/2:Nxi/2-1)*dxi/lambda; % fftshifted frequency-plane coordinates.  dimensions = 1xN
            FX2 = ones(Nxi,1)*fx.^2; % squared, shifted frequency-plane coordinates. dimensions = NxN
            H = exp(-1i*pi*lambda*dz*(FX2+FX2.')); % transfer function
            
            %--Pupil Plane Coordinates
            if( mod(M,2)==1 )
                xs = ( -(M-1)/2:(M-1)/2 ).'*dx;
            elseif( (mod(M,2)==0) && strcmpi(centering,'interpixel') )
                xs = ( -(M-1)/2:(M-1)/2 ).'*dx;
            else
                xs = (-M/2:(M/2-1)).'*dx;
            end
            
            %--Focal Plane Coordinates
            if(  (mod(Nxi,2)==1) ) %--Odd-sized array
                xis = ( -(Nxi-1)/2:(Nxi-1)/2 )*dxi;
            elseif(strcmpi(centering,'interpixel'))%--Even-sized array, interpixel centered
                xis = ( -(Nxi-1)/2:(Nxi-1)/2 )*dxi;
            else%--Even-sized array, pixel centered
                xis = ( -Nxi/2:(Nxi/2-1) )*dxi;
            end

            rect_mat_post  = (exp(-2*pi*1i*(xs*xis)/(lambda*f))); %--Just the transpose of rect_mat_pre for square input and output matrices with same samplings in x and y.
            U1 = 1/(lambda*f)*(rect_mat_post.'*Ein*rect_mat_post);
            Eout = (dx*dx)*(dxi*dxi)/(lambda*f)*conj(rect_mat_post)*(H.*U1)*rect_mat_post';

    end
end %--END OF FUNCTION