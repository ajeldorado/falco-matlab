function phz = falco_zwfs_reconstructor(I0,IZ, mask, b, theta, type)
% phz = falco_zwfs_reconstructor(I0,IZ, mask, b, theta, type)
%   Function to estimate the phase from Zernike WFS data
%   
%   Inputs:
%       I0 - image without Zernike mask
%       IZ - image with Zernike mask aligned 
%       mask - pupil support mask 
%       b - Reference wave
%       theta - mask phase shift, typically 2*pi*depth*(n-1)/lambda;
%       type - string - 'linear', 'ndiaye', or 'full' reconstructors 
%
%   Output:
%       phz - A 2D array with the estimated phase 

    A = sqrt(I0);
	b = mean(A(mask))*b; % Scale reference wave to match incident energy
    
	if( strcmpi(type,'linear') || strcmpi(type,'l') ) % Linear approximation 
        
        denom = 2.*b.*A.*sin(theta);
        term2 = A.^2 + 2*b.^2 - 2*b.*A + 2*b.*(A-b).*cos(theta);
        phz = (IZ - term2)./(denom+1e-30);
        
    elseif( strcmpi(type,'ndiaye') || strcmpi(type,'n') ) % Second order approximation

        a1 = 2*b.*A.*(-2*b.^2+IZ+3*b.*A-A.^2+b.*(2*b-A)*cos(theta))*sin(theta/2)^2;
        a2 = b.*A.*sin(theta);
        numer = sqrt(a1) + a2;
        denom = b.*A*(cos(theta)-1);
        phz = -1*real(numer./(denom+1e-30)); 
        
    elseif( strcmpi(type,'full') || strcmpi(type,'f') ) % Full analytical solution 
        
        a0 = 4*b.^2.*(IZ+A.^2)-6*b.^4-(IZ-A.^2).^2-4*b.^2.*(IZ+A.^2-2*b.^2)*cos(theta)-2*b.^4*cos(2*theta);
        a1 = b.^2.*A.^2.*sin(theta)^2.*a0;
        a2 = b.*A.*(A.^2-IZ).*(cos(theta)-1)-2*b.^3.*A.*(cos(theta)-1)^2;
        numer = a2 + sqrt(a1); % plus or minus are valid solutions 
        denom = 4.*b.^2.*A.^2.*(cos(theta)-1);
        phz = real(acos(numer./(denom+1e-30))); % plus or minus are valid solutions 
        
    else
        error('Reconstructor undefined');
    end

	phz = phz - mean(phz(mask));

end