function phz = falco_zwfs_reconstructor(IC0,IC, mask, b, type,varargin)
%phz = falco_zwfs_reconstructor(IC0,IC, mask, b, type, (modelparams))
%   IC0 - image without Zernike mask
%   IC - iamge with Zernike mask aligned 
%   mask - pupil support mask 
%   b - Reference wave
%   type - string - 'wallace' or 'ndiaye' reconstructors 
%   Need to pass mp when using ndiaye

    A02 = mean(IC0(mask));
    P = sqrt(IC0/A02);

    if(strcmp(type,'wallace')||strcmp(type,'Wallace')||strcmp(type,'w')||strcmp(type,'W'))
        numer = IC - A02*(P.^2 + 2*b.^2);
        denom = 2*sqrt(2)*A02*P.*b;
        phz = real(pi/4 + asin(numer./denom));
    elseif(strcmp(type,'NDiaye')||strcmp(type,'Ndiaye')||strcmp(type,'ndiaye')||strcmp(type,'n')||strcmp(type,'N'))
        mp = varargin{1};
        theta = 2*pi*mp.F3.t*(mp.F3.n(mp.lambda0)-1)/mp.lambda0;
        a1 = 2*b.*P.*(-2*b.^2+IC+3*b.*P-P.^2+b.*(2*b-P)*cos(theta))*sin(theta/2)^2;
        a2 = b.*P.*sin(theta);
        numer = sqrt(a1) + a2;
        denom = b.*P.*(cos(theta)-1);
        phz = real(numer./(denom+1e-30)); 
    else
        disp('Reconstructor undefined');
    end

	phz = phz - mean(phz(mask));

end

