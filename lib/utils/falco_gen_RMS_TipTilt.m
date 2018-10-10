% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--OVERVIEW:
% Function to generate the weights and (x,y) coordinates of the intensity
% for tip/tilt from stellar angular size and pointing jitter.
%
%--INPUTS:
% TTrms:    RMS tip/tilt value, in milliarcseconds
% Dstar:    stellar angular diameter, in milliarcseconds
% Dtel:     diameter of the telescope, in meters (used for conversion from mas to lambda0/D)
% lambda:    center wavelength, in meters (used for conversion from mas to lambda0/D)
% 'Nacross',value:  (optional, needs keyword) Number of points across the 2-D weight array (want as odd number)
% 'Rfac',value:     Zero out weights outside this cutoff radius
% 'upsampFac',value: Upsampling factor for computing the resolution at
%     which to convolve the stellar disk and the Gaussian tip/tilt jitter profile.
%
%--OUTPUTS:
% xsTT:  Vector of x-coordinates for the jitter weights, in lambda0/D
% ysTT:  Vector of y-coordinates for the jitter weights, in lambda0/D
% wsTT:  Vector of jitter weights (the sum of the weights is one)
%
%--REVISION HISTORY: 
% Created by A.J. Riggs on 2017-11-01.
%
%--EXAMPLE FUNCTION CALL:
% [rmsTT_x_vec,rmsTT_y_vec,rmsTT_w_vec] = falco_gen_RMS_TipTilt(TTrms,Dstar,Dtel,lambda0,'Nacross',7,'Rfac',3,'upsampFac',256)


function [xsTT,ysTT,wsTT] = falco_gen_RMS_TipTilt(TTrms,Dstar,Dtel,lambda,varargin)

if(TTrms==0 && Dstar==0)
    xsTT = 0;
    ysTT = 0;
    wsTT = 1;
else

% Set default values of input parameters
%  7,3 is good for real evaluation (29 points), and 5,2 is good for quicker real-time evaluation (13 points)
Nacross = 7; %--Number of points across the 2-D weight array (want as odd number)
Rfac = 3;% 3; %--Zero out weights outside this cutoff radius.    
upsampFac = 256;

flagBigD = TTrms<0.25*Dstar;  % Case 1: Stellar diameter dominates

icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
icav = icav + 1;
    switch lower(varargin{icav})
        case {'nacross'}
            icav = icav + 1;
            Nacross   = varargin{icav};  % aperture center to beam center X
        case {'rfac'}
            icav = icav + 1;
            Rfac   = varargin{icav};  % aperture center to beam center Y
        case {'upsampfac'}
            icav = icav + 1;
            upsampFac = varargin{icav};
        otherwise
            error('falco_gen_RMS_TipTilt: Unknown keyword: %s\n', varargin{icav});
    end
end

  
mas2lam0D = 1/(lambda/Dtel*180/pi*3600*1000);
NacrossHD = Nacross*upsampFac;


%--Coordinates (high-res for the convolution, and low-res for the output)
%  Coordinates are pixel-centered on an odd-sized array.

if(flagBigD) % Case 1: Stellar diameter dominates
    if(mod(Nacross,2)==1)
        dx = Dstar/(Nacross-1);
    else
        dx = Dstar/Nacross;
    end
else
    dx = TTrms; % Case 2: RMS pointing jitter dominates
end


xs = (-(Nacross-1)/2:1:(Nacross-1)/2)*dx;
xsHD = (-(NacrossHD-1)/2:1:(NacrossHD-1)/2)*(dx/upsampFac);

[JXS,JYS] = meshgrid(xs);
RS = sqrt(JXS.^2 + JYS.^2);

[JXShd,JYShd] = meshgrid(xsHD);
RShd = sqrt(JXShd.^2 + JYShd.^2);

%--Define the stellar disk at high resolution
stellar_disk_hd = 0*RShd;
stellar_disk_hd(RShd<=Dstar/2) = 1;
stellar_disk = 0*RS;
for kk=1:Nacross
    for ll=1:Nacross
        stellar_disk(kk,ll) = sum(sum(stellar_disk_hd((kk-1)*upsampFac+1:kk*upsampFac,(ll-1)*upsampFac+1:ll*upsampFac)))/upsampFac^2;
    end
end

% figure(91); imagesc(xsHD,xsHD,stellar_disk_hd); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% % figure(92); imagesc(xs,xs,stellar_disk); axis xy equal tight; colorbar;set(gca,'Fontsize',20);


mask = ones(Nacross);
if(flagBigD) 
    mask(RS > 1.2*Dstar/2) = 0; % Case 1: Stellar diameter dominates
else
    mask(RS > Rfac*TTrms) = 0; % Case 2: RMS pointing jitter dominates
end


%     figure(301); imagesc(xs,xs,mask);
j_ele = find(mask==1);
Njitter = sum(sum(mask));

jw_mat_hd = exp(-1/2*(RShd/TTrms).^2); %--Matrix of jitter weights (jw)
% figure(93); imagesc(xsHD,xsHD,jw_mat_hd); axis xy equal tight; colorbar; set(gca,'Fontsize',20);

%--Convolve Gaussian for the jitter with the disk for stellar angular size. (Both at high-res.)
if(Dstar==0)
    jw_mat_hd_sas = jw_mat_hd;
elseif(TTrms==0)
    jw_mat_hd_sas = stellar_disk_hd;
else
    jw_mat_hd_sas = ifftshift(ifft2( fft2(fftshift(jw_mat_hd)).*fft2(fftshift(stellar_disk_hd)) ));
end

% figure(94); imagesc(xsHD,xsHD,jw_mat_hd_sas); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
%figure(95); imagesc(xsHD,xsHD,jw_mat_hd_sas-jw_mat_hd); axis xy equal tight; colorbar; set(gca,'Fontsize',20);

%--Bin down the T/T weights to the sampling desired.
jw_mat = 0*RS;
for kk=1:Nacross
    for ll=1:Nacross
        jw_mat(kk,ll) = sum(sum(jw_mat_hd_sas((kk-1)*upsampFac+1:kk*upsampFac,(ll-1)*upsampFac+1:ll*upsampFac)))/upsampFac^2;
    end
end

jw_mat = jw_mat/max(jw_mat(:));
%     figure(302); imagesc(xs,xs,jw);

jw_mat = jw_mat.*mask;
jw_mat = jw_mat/sum(sum(jw_mat));

wsTT = jw_mat(j_ele).'; %--Vector of jitter weights
xsTT = mas2lam0D*JXS(j_ele).'; %--Vector of x-coordinates for the jitter weights, in lambda0/D
ysTT = mas2lam0D*JYS(j_ele).'; %--Vector of y-coordinates for the jitter weights, in lambda0/D
% figure(303); imagesc(xs,xs,jw_mat); colorbar; 
 
end

end %--END OF FUNCTION

