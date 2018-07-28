% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = propcustom_PTP_for_GPU(Ein,FR2shift,lambda,dz)
%
% This is a Fresnel propagation function. It uses the convolution, transfer function form of the Fresnel integral.  
%
% This version requires a matrix with an even, equal number of points in each dimension.
%
%
% REQUIRED INPUTS: 
%  Ein = input E-field as 2-D square array
%  N = number of points along both dimensions of Ein
%  dx = linear resolution at each plane (meters/pixel)
%  lambda = wavelength of light (meters)
%  dz = propagation distance (meters)
%
% OUTPUTS:
%  Eout: output E-field as square array
%
% Modified on 2018-07-27 by A.J. Riggs not to use meshgrid. Loads FR2shift.
% Modified on 2018-04-26 by A.J. Riggs to correct the sampling requirement.
% Written by A.J. on April 30, 2013. 


% Within the GPU:
function Eout = propcustom_PTP_for_GPU(Ein,N,dx,lambda,dz)


    %--Make coordinate vectors and maps in the frequency plane
    fx = [0:(N/2-1),(-N/2:-1)]/(N*dx);    % fftshifted frequency-plane coordinates.  dimensions = 1xN
    FX2 = ones(N,1)*fx.^2; % squared, shifted frequency-plane coordinates. dimensions = NxN

    H =  exp(-1i*pi*lambda*dz*(FX2 + FX2.')); % transfer function. can use if exp() function exists in CUDA.
    U1=fft2(fftshift(Ein));      % shift, then FFT source field
    Eout=ifftshift(ifft2(H.*U1));    % element-wise multiply, inverse FFT, and then unshift    

end


% %--Make coordinate vectors and maps in the frequency plane
% fx = ( -N/2:1:(N/2-1) )/(N*dx);    % frequency-plane coordinates.  dimensions = 1xN
% fx2 = fx.^2; % squared frequency-plane coordinates. dimensions = 1xN
% FX2 = ones(N,1)*fx2;
% FR2shift = fftshift(FX2 + FX2.');
% figure(1); imagesc(FX2); axis xy equal tight; colorbar;
% figure(2); imagesc(FR2shift); axis xy equal tight; colorbar;









% % REQUIRED INPUTS: 
% %  Ein = input E-field as 2-D square array
% %  FR2shift = 2-D quadratic phase term in convolution (1/meters)
% %  lambda = wavelength of light (meters)
% %  dz = propagation distance (meters)
% %
% % OUTPUTS:
% %  Eout: output E-field as square array
% %
% % Modified on 2018-07-27 by A.J. Riggs not to use meshgrid. Loads FR2shift.
% % Modified on 2018-04-26 by A.J. Riggs to correct the sampling requirement.
% % Written by A.J. on April 30, 2013. 
% 
% 
% % Within the GPU:
% function Eout = propcustom_PTP_for_GPU(Ein,FR2shift,lambda,dz)
%     
%     H =  exp(-1i*pi*lambda*dz*FR2shift); % transfer function. can use if exp() function exists in CUDA.
%     U1=fft2(fftshift(Ein));      % shift, then FFT source field
%     Eout=ifftshift(ifft2(H.*U1));    % element-wise multiply, inverse FFT, and then unshift    
% 
% end