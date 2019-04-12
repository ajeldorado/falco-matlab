% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function[Zern_coeff,Amat] = propcustom_getZernCoefAndMat_masked(Phase_in,minZerN,maxZerN,mask)

% fprintf('Initializing Zernikes...');
M = size(Phase_in,1);
xpup = -1:(2/(M-1)):1; %spacing is just constant
ypup = -1:(2/(M-1)):1;
k = 1;

for i = 1:length(xpup)
    for j = 1:length(ypup)
        if((xpup(i)^2 + ypup(j)^2)<= 1)
            [th(k), r(k)] = cart2pol(xpup(i), ypup(j));
            xVec(k) = i;
            yVec(k) = j;
            k = k+1;
        end
    end
end
indx = sub2ind([M M], xVec, yVec);

% number of Zernike modes
sz = .5*(maxZerN+1)*(maxZerN+2) - .5*(minZerN+1)*(minZerN);
A = zeros(M*M, sz);
oddeven = zeros(sz, 1);
% Evaluate and store Zernike polynomials
k = 1;
for n = minZerN:maxZerN
    for m = -n:2:n
        ZP = zeros(M,M);
        ZP(indx) = propcustom_zernfun(n, m, r, th, 'norm');
        A(:, k) = reshape(ZP.*mask,M*M,1);
        
        if rem(m, 2) ~= 0 
            oddeven(k) = 1;
        end
        k = k + 1;
    end
end

B = reshape(Phase_in,M*M,1);

Zern_coeff = pinv(A)*B;
Amat = reshape(A,[M,M,size(A,2)]);
end