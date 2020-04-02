function phz = falco_zwfs_diff_reconstructor(I0, IZ, Iref, mask, b, theta, type)
% phz = falco_zwfs_diff_reconstructor(I0, IZ, Iref, mask, b, theta, type)
%   Function to estimate the phase difference from Zernike WFS measurements
%   
%   Inputs:
%       I0 - image with Zernike mask misaligned 
%       IZ - image with Zernike mask aligned 
%       Iref - Reference image (also with Zernike mask aligned)
%       mask - pupil support mask 
%       b - Reference wave
%       theta - mask phase shift, typically 2*pi*mp.F3.t*(mp.F3.n(mp.lambda0)-1)/mp.lambda0;
%       type - string - 'linear' reconstructor 
%
%   Output:
%       phz - A 2D array with the estimated phase 

    P = sqrt(I0);
    b = mean(P(mask))*b; % Scale reference wave to match incident energy

	if( strcmpi(type,'linear') || strcmpi(type,'l') )
        
        denom = 2*b.*P.*sin(theta);
        phz = (IZ - Iref)./denom;
        
    else
        error('Reconstructor undefined');
    end

	phz = phz - mean(phz(mask));

end