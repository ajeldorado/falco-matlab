% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------


function Aout = falco_bin_downsample(Ain,dsfac)
%     Downsample an array by binning.
% 
%     Parameters
%     ----------
%     Ain : 2-D array
%         The matrix to be downsampled
%     dsfac : int
%         Downsampling factor for the matrix
% 
%     Returns
%     --------
%     Aout : 2-D array
%         Downsampled array

    % Error checks on inputs
    if(length(size(Ain))~=2); error('Ain must be a matrix.'); end
    if(mod(dsfac,1)~=0); error('Argument dsfac must be an integer.'); end
    
    % Array Sizes
    Nx0 = size(Ain,2);
    Ny0 = size(Ain,1);
    if(mod(Nx0,dsfac)~=0 || mod(Ny0,dsfac)~=0  )
        error('The size of Ain must be divisible by dsfac.');
    end
    Nx1 = Nx0/dsfac;
    Ny1 = Ny0/dsfac;
    
    % Bin and average values from the high-res array into the low-res array
    Aout = zeros(Ny1,Nx1);
    for ix = 1:Nx1
        for iy = 1:Ny1
            Aout(iy,ix) = sum(sum(Ain(dsfac*(iy-1)+1:dsfac*iy,dsfac*(ix-1)+1:dsfac*ix)))/dsfac^2;
        end
    end

end


