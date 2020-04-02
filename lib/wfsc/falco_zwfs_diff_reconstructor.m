function phz = falco_zwfs_diff_reconstructor(I0, IZ, Iref, mask, b, theta, varargin)
% phz = falco_zwfs_diff_reconstructor(I0, IZ, Iref, mask, b, theta, <phz_ref>)
%   Function to estimate the phase difference from Zernike WFS measurements
%   Assumes that the difference in phase between IZ and Iref is small. 
%   Also allows for a correction based on an estimate of the phase 
%   corresponding to Iref. 
%   
%   Inputs:
%       I0 - image with Zernike mask misaligned 
%       IZ - image with Zernike mask aligned 
%       Iref - Reference image (also with Zernike mask aligned)
%       mask - pupil support mask 
%       b - Reference wave
%       theta - mask phase shift, typically 2*pi*mp.F3.t*(mp.F3.n(mp.lambda0)-1)/mp.lambda0;
%       varargin{1} - phz_ref - The reference absolute phase estimate 
%
%       If either the phase difference is substantial,
%       use the full reconstructor in falco_zwfs_reconstructor
%
%   Output:
%       phz - A 2D array with the estimated phase 

    A = sqrt(I0);
    b = mean(A(mask))*b; % Scale reference wave to match incident energy
    
    if(nargin==1)
        phz_ref = varargin{1};
        beta = sin(theta)*cos(phz_ref) + (1-cos(theta))*sin(phz_ref);
    else
        beta = sin(theta);
    end
    
    denom = 2*b.*A.*beta;
    phz = (IZ - Iref)./denom;% first order phase difference approximation 

	phz = phz - mean(phz(mask));% subtract of constant offset 

end