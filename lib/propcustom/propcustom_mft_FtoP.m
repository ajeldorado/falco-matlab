% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function [Epup] = propcustom_mft_FtoP(Efoc,fl,lambda,dxi,deta,dx,N,varargin)
%--This function propagates from the focus to a pupil via a Matrix Fourier Transform (MFT).
%
%--PROPERTIES
% -Unconventional normalization used so that energy is conserved when
% squaring and summing the E-field directly (without multiplying by the
% dx*dy in Parseval's Theorem).
%
%--INPUTS
% Efoc = electric field at input (focal) plane
% fl = focal length of focusing optic in meters
% lambda = wavelength in meters
% dxi = pixel width in focal plane (CCD pixel size) in meters
% eta = pixel height in focal plane (CCD pixel size) in meters
% dx = resolution of a pixel in the pupil plane in meters
% N = number of points across the pupil plane
% 'INTERPIXEL' = optional flag to change the centering of the array to be between pixels
%
% - Modified on 2019-04-05 by A.J. Riggs to remove the 1/1i term from each FT.
%
%--------------------------------------------------------------------------

function [Epup] = propcustom_mft_FtoP(Efoc,fl,lambda,dxi,deta,dx,N,varargin)

% Set default value(s) of the optional input parameter(s)
centering = 'pixel';

%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'interpixel'}
            centering = 'interpixel'; % For even arrays, beam center is in between pixels.
        case {'pixel'}
            centering = 'pixel'; % For even arrays, beam center is in between pixels.    
        otherwise
            error('propcustom_mft_PtoF: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

%--Focal Plane Coordinates
[Neta,Nxi] = size(Efoc);
if(  (mod(Nxi,2)==1) ) %--Odd-sized array
    xis = ( -(Nxi-1)/2:(Nxi-1)/2 ).'*dxi;
elseif(strcmpi(centering,'interpixel'))%--Even-sized array, interpixel centered
    xis = ( -(Nxi-1)/2:(Nxi-1)/2 ).'*dxi;
else%--Even-sized array, pixel centered
    xis = ( -Nxi/2:(Nxi/2-1) ).'*dxi;
end

if(  (mod(Neta,2)==1) ) %--Odd-sized array
    etas = ( -(Neta-1)/2:(Neta-1)/2 )*deta;
elseif(strcmpi(centering,'interpixel'))%--Even-sized array, interpixel centered
    etas = ( -(Neta-1)/2:(Neta-1)/2 )*deta;
else%--Even-sized array, pixel centered
    etas = ( -Neta/2:(Neta/2-1) )*deta;
end

%--Pupil Plane Coordinates
if( mod(N,2)==1 )
    %dx = L2/(N-1);
    xs = ( -(N-1)/2:(N-1)/2 )*dx;
elseif( (mod(N,2)==0) && strcmpi(centering,'interpixel') )
    %dx = L2/N;
    xs = ( -(N-1)/2:(N-1)/2 )*dx;
else
    %dx = L2/N;
    xs = (-N/2:(N/2-1))*dx;
end
ys = xs.';
dy = dx;

%--Matrix Fourier Transform (MFT)
rect_mat_pre = (exp(-2*pi*1i*(ys*etas)/(lambda*fl)));
rect_mat_post  = (exp(-2*pi*1i*(xis*xs)/(lambda*fl)));
Epup = sqrt(dx*dy)*sqrt(dxi*deta)/(1*lambda*fl)*(rect_mat_pre*Efoc*rect_mat_post); 

end %--END OF FUNCTION