% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function [Efoc] = propcustom_mft_PtoF(Epup, f,lambda,dx,dxi,Nxi,deta,Neta,varargin)
%
%--This function propagates from a pupil to a focus via a Matrix Fourier Transform (MFT).
%
%--PROPERTIES
% -Unconventional normalization used so that energy is conserved when
% squaring and summing the E-field directly (without multiplying by the
% dxi*deta in Parseval's Theorem).
%
%--INPUTS
% Epup = electric field at input (pupil) plane
% fl = focal length of focusing optic in meters
% lambda = wavelength in meters
% dx = resolution of a pixel in the pupil plane in meters
% dxi = pixel width in focal plane (CCD pixel size) in meters
% Nxi = total number of points horizontally across at focal plane (second plane)
% eta = pixel height in focal plane (CCD pixel size) in meters
% Neta = total number of points vertically across at focal plane (second plane)
%
%--OPTIONAL INPUTS (via keywords)
% 'INTERPIXEL' = flag to change the centering of the array to be between pixels
% 'xpc' = change the x-center of the pupil plane to the value after this flag
% 'ypc' = change the y-center of the pupil plane to the value after this flag
% 'xfc' = change the x-center of the focal plane to the value after this flag
% 'yfc' = change the y-center of the focal plane to the value after this flag
%--------------------------------------------------------------------------

function [Efoc] = propcustom_mft_PtoF(Epup, f,lambda,dx,dxi,Nxi,deta,Neta,varargin)

% Set default value(s) of the optional input parameter(s)
centering = 'pixel';
xpc = 0;
ypc = 0;
xfc = 0;
yfc = 0;

%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'interpixel'}
            centering = 'interpixel'; % For even arrays, beam center is in between pixels.
        case {'pixel'}
            centering = 'pixel'; % For even arrays, beam center is in between pixels.    
        case {'xpc'} 
            icav = icav + 1;
            xpc = varargin{icav};  % x-center location in pupil plane [meters]
        case {'ypc'} 
            icav = icav + 1;
            ypc = varargin{icav};  % y-center location in pupil plane [meters]
        case {'xfc'} 
            icav = icav + 1;
            xfc = varargin{icav};  % x-center location in focal plane [meters]
        case {'yfc'} 
            icav = icav + 1;
            yfc = varargin{icav};  % y-center location in focal plane [meters]
            
        otherwise
            error('propcustom_mft_PtoF: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end



%--Pupil Plane Coordinates
[M,N] = size(Epup);
if M~=N % Just use square inputs for the time being. Can change later
    disp('Error: input matrix is not square');
end
if( mod(M,2)==1 )
    %dx = L1/(M-1);
    xs = ( -(M-1)/2:(M-1)/2 ).'*dx;
elseif( (mod(M,2)==0) && strcmpi(centering,'interpixel') )
    %dx = L1/M;
    xs = ( -(M-1)/2:(M-1)/2 ).'*dx;
else
    %dx = L1/M;
    xs = (-M/2:(M/2-1)).'*dx;
end
ys = xs.';
dy = dx;

%--Translate the pupil plane coordinates
xs = xs - xpc;
ys = ys - ypc;


%--Focal Plane Coordinates
if(  (mod(Nxi,2)==1) ) %--Odd-sized array
    xis = ( -(Nxi-1)/2:(Nxi-1)/2 )*dxi;
elseif(strcmpi(centering,'interpixel'))%--Even-sized array, interpixel centered
    xis = ( -(Nxi-1)/2:(Nxi-1)/2 )*dxi;
else%--Even-sized array, pixel centered
    xis = ( -Nxi/2:(Nxi/2-1) )*dxi;
end

if(  (mod(Neta,2)==1) ) %--Odd-sized array
    etas = ( -(Neta-1)/2:(Neta-1)/2 ).'*deta;
elseif(strcmpi(centering,'interpixel'))%--Even-sized array, interpixel centered
    etas = ( -(Neta-1)/2:(Neta-1)/2 ).'*deta;
else%--Even-sized array, pixel centered
    etas = ( -Neta/2:(Neta/2-1) ).'*deta;
end

%--Translate the focal plane coordinates
xis = xis - xfc;
etas = etas - yfc;


%--Matrix Fourier Transform (MFT)
rect_mat_pre = (exp(-2*pi*1i*(etas*ys)/(lambda*f)));
rect_mat_post  = (exp(-2*pi*1i*(xs*xis)/(lambda*f)));
Efoc = sqrt(dx*dy)*sqrt(dxi*deta)/(1i*lambda*f)*(rect_mat_pre*Epup*rect_mat_post);

end %--END OF FUNCTION
