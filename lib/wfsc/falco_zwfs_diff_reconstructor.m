function phz = falco_zwfs_diff_reconstructor(I0, IZ, Iref, mask, b, theta, type, varargin)
% phz = falco_zwfs_diff_reconstructor(I0, IZ, Iref, mask, b, theta, type)
%   Function to estimate the phase difference from Zernike WFS measurements
%   Assumes that the difference in phase between I0 and Iref is small. 
%   Allows for a correction based on an estimate of the phase corresponding
%   to Iref. 
%   
%   Inputs:
%       I0 - image with Zernike mask misaligned 
%       IZ - image with Zernike mask aligned 
%       Iref - Reference image (also with Zernike mask aligned)
%       mask - pupil support mask 
%       b - Reference wave
%       theta - mask phase shift, typically 2*pi*mp.F3.t*(mp.F3.n(mp.lambda0)-1)/mp.lambda0;
%       type - string - 'linear' for linear reconstructor, 'first' for first order correction 
%                       anything beyond linear requires the reference phase estimate 
%                       phz_ref to be passed as first varargin. 
%
%       Use 'linear' when you have no clue what the phase is for the
%       reference image, but you know it is small (<< 1 radian). 
%       Introducing an estimate of phz_ref can reduce this error.
%
%       If either the reference phase or phase difference are substantial,
%       use the full reconstructor in falco_zwfs_reconstructor
%
%   Output:
%       phz - A 2D array with the estimated phase 

    P = sqrt(I0);
    b = mean(P(mask))*b; % Scale reference wave to match incident energy
    
    denom = 2*b.*P.*sin(theta);
    phz = (IZ - Iref)./denom;% first order phase difference approximation 

	if( strcmpi(type,'linear') || strcmpi(type,'l') )
        
        % we're done here
        correction = 1; 
        
	elseif( strcmpi(type,'first') || strcmpi(type,'firstorder') )
        
        phz_ref = varargin{1}; 
        
        a1 = (1-cos(theta))/sin(theta);
        
        correction = 1 - a1*phz_ref; 
        
	elseif( strcmpi(type,'second') || strcmpi(type,'secondorder') )
        
        phz_ref = varargin{1}; 
        
        a1 = (1-cos(theta))/sin(theta);
        a2 = (1/cos(theta/2)^2 - 1/2);
        
        correction = 1 - a1*phz_ref + a2*phz_ref.^2; 
        
	else
        error('Reconstructor undefined');
    end
    
    phz = phz.*correction; 

	phz = phz - mean(phz(mask));

end