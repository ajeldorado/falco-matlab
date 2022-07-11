% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
function [phz, phz_recon] = removeZernikes(phz, N, M, ROI)
%[phz, phz_recon] = removeZernikes(phz, N, M, ROI)
%    Removes Zernike terms from phase.
%
%    Inputs: 
%        phz - the input phase. 
%        N,M - list of n,m indices (i.e. Z_n^m) to remove.
%        ROI - Region of interest.
%    Returns: 
%        phz - phz with orders removed 
%        phz_recon - the reconstructed phase 

    [rows, cols, dims] = size(phz); 
    
    [X, Y] = meshgrid(-cols/2:cols/2-1, -rows/2:rows/2-1);
    [THETA, RHO] = cart2pol(X, Y);
    Z = zernfun(N, M, RHO(ROI)/max(max(RHO(ROI))), THETA(ROI));
    
    phz_recon = zeros(size(phz));
    for slice_num = 1:dims
       
        im = phz(:, :, slice_num);
        a = Z\im(ROI);
        reconstructed = NaN(size(im));
        reconstructed(ROI) = Z*a; 
        phz_recon(:, :, slice_num) = reconstructed;
        
    end
    
    phz = phz - phz_recon;
    phz(isnan(phz)) = 0;
    phz_recon(isnan(phz_recon)) = 0;

end

