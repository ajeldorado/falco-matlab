% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a datacube of Zernike modes with 3-fold symmetry.
%
% Modified to be a function on 2018-11-20 by A.J. Riggs
% Created on 2016-11-21 by A.J. Riggs.

function ZmapCube = falco_3fold_symmetry_Zernikes(Nbeam,maxRadialOrder,centering,varargin)


% Set default values of input parameters
  symmaxis = 'none';            %--axis of mirror symmetry ('x', 'y', or default of 'none')

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'symmaxis'}
        icav = icav + 1;
        symmaxis   = varargin{icav};  % input array center X (pixels)
      otherwise
        error('falco_3fold_symmetry_Zernikes: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  
% %% DEBUGGING ONLY: Hard-Coded values for Debugging as a Standalone Script
% clear all; 
% 
% Nbeam = 100;
% symmaxis = 'y';
% centering = 'pixel';
% maxRadialOrder = 10;%22;%20;%30;%10;%5;%4;%5;%15;
% 
% addpath('~/Repos/falco-matlab/lib/external'); 
% addpath('~/Repos/falco-matlab/lib/Zernike/');
% addpath('~/Documents/MATLAB/PROPER/');
% %%


%--LEAVE THIS ALONE
method = 'zernfun'; %'PROPER';%'zernfun';

Q = floor(maxRadialOrder/3);  % Number of 3-fold mode types (3,6,9, etc, not including 0).
%--Calculate number of Zernike modes kept:
Nmodes = floor(maxRadialOrder/2)+1; % 0th order modes
for q=1:Q
    Nmodes = Nmodes + 2*floor((maxRadialOrder-3*q+2)/2);
end
fprintf('Max Zernike order = %d.\t\t Number of 3-fold symmetric modes = %d\n',maxRadialOrder,Nmodes);

%--Coordinates normalized to the beam radius (not diameter)
switch centering
    case 'interpixel'
        Narray = ceil_even(Nbeam);
        xs = (-(Narray-1)/2:(Narray-1)/2)/Nbeam*2;
    case 'pixel'
        Narray = ceil_even(Nbeam+1);
%         xs = (-(Narray-1)/2:(Narray-1)/2)/Nbeam*2;
        xs = (-Narray/2:(Narray/2-1))/Nbeam*2;
end
[XS,YS] = meshgrid(xs);
RS = sqrt(XS.^2+YS.^2); 
THETAS = atan2(YS,XS);
mask = RS<=1;%
% figure; imagesc(mask); colorbar;


%%
switch method
%% Method using zernfun.m to generate Zernike maps. Limitation: can't set maxRadialOrder>22.

    case 'zernfun'

Phase_in = zeros(Narray);


% fprintf('Initializing Zernikes...');
M = size(Phase_in,1);

% switch centering
%     case 'interpixel'
%         xpup = -1:(2/(M-1)):1; %spacing is just constant
%         ypup = -1:(2/(M-1)):1;
%     case 'pixel'
%         xpup = (-M/2:M/2-1)/(M/2);%-1:(2/(M-1)):1; %spacing is just constant
%         ypup = (-M/2:M/2-1)/(M/2);%-1:(2/(M-1)):1;   
% end

xpup = xs;
ypup = xs;

kk = 1;
for ii = 1:Narray %length(xpup)
    for jj = 1:Narray %length(ypup)
        if((xpup(ii)^2 + ypup(jj)^2)<= 1)
            %[th(kk), r(kk)] = cart2pol(xpup(ii), ypup(jj));
            th(kk) = THETAS(jj,ii);
            r(kk) = RS(jj,ii);
            xVec(kk) = ii;
            yVec(kk) = jj;
            kk = kk+1;
        end
    end
end
indx = sub2ind([M M], xVec, yVec);
% fprintf('done.\n');
% whos indx

% Evaluate and store Zernike polynomials

fprintf('Making datacube of Zernikes...');
Nmodes = Nmodes;%.5*(maxRadialOrder+1)*(maxRadialOrder+2) - .5*(minZerN+1)*(minZerN); % number of Zernike modes
% sz = .5*(maxRadialOrder+1)*(maxRadialOrder+2) - .5*(minZerN+1)*(minZerN); % number of Zernike modes
ZmapVec = zeros(M*M, Nmodes);
% oddeven = zeros(sz, 1);



% Circularly symmetric modes only (angular frequency m=0 )
kk = 1; % A vector index
for n = 0:2:maxRadialOrder
    m = 0;
    ZP = zeros(M,M);
    ZP(indx) = zernfun(n, m, r, th, 'norm');
    ZmapVec(:, kk) = reshape(ZP,M*M,1);
% 
%     if rem(m, 2) ~= 0 
%         oddeven(k) = 1;
%     end
    kk = kk + 1;
    
end
%


% 3x-order Modes: 
minpos = [-1, 1];
for q=1:Q
    for n = 3*q:2:maxRadialOrder
    %     fprintf('%d ',n);
        for pm = 1:2 %m = -n:2:n
            m = 3*q*minpos(pm);
%             fprintf('n,m = \t%d,%d\n',n,m);
            ZP = zeros(M,M);
            ZP(indx) = zernfun(n, m, r, th, 'norm');
            ZmapVec(:, kk) = reshape(ZP,M*M,1);
            
%             Atemp = reshape(A(:,k),[M,M]); %--for debugging
%             figure(20); imagesc(Atemp); axis xy equal tight;
%             title(sprintf('Mode %d',k),'Fontsize',20);
%             pause(0.5);
            
%             if rem(m, 2) ~= 0 
%                 oddeven(k) = 1;
%             end
            kk = kk + 1;
        end
    end
end


% fprintf('Making fit of wavefront to Zernikes...');
% B = reshape(Phase_in,M*M,1);
% Zern_coeff = pinv(ZmapVec)*B;%(A'*A)\(A'*B);

ZmapCubeInit = reshape(ZmapVec,[M,M,Nmodes]);
fprintf('done.\t');

% for ii=1:Nmodes %floor(maxRadialOrder/2)+1
%     figure(20); imagesc(ZmapCubeInit(:,:,ii)); axis xy equal tight;
%     title(sprintf('Mode %d',ii),'Fontsize',20);
%     pause(1/20);
% end


%% Method using PROPER to generate Zernikes. Limitation: can't set maxRadialOrder>22.

    case 'PROPER'
        
indsZnoll = zeros(Nmodes,1); % []; %--Initialize the Noll indices of Zernikes to compute values for
zz = 1;

% Circularly symmetric modes only (angular frequency m=0 )
for n = 0:2:maxRadialOrder
    m = 0;
    
    %--Compute the Zernike Noll index
    if( (mod(n,4)==0) || (mod(n,4)==1) )
        if(m>0)
            val = 0;
        else
            val = 1;
        end
    else % elseif( (mod(n,4)==2) || (mod(n,4)==3) )
        if(m<0)
            val = 0;
        else
            val = 1;
        end
    end
    Znoll = n*(n+1)/2+abs(m)+val; %--index of Zernike Noll indices used

    indsZnoll(zz) = Znoll;
%     indsZnoll = [indsZnoll; Znoll];
    zz = zz+1;
end 


% 3-fold Symmetry Zernike Modes: 
minpos = [-1, 1];
for q=1:Q
    for n = 3*q:2:maxRadialOrder
    %     fprintf('%d ',n);
        for pm = 1:2 %m = -n:2:n
            m = 3*q*minpos(pm);
%             fprintf('n,m = \t%d,%d\n',n,m);

            %--Compute the Zernike Noll index
            if( (mod(n,4)==0) || (mod(n,4)==1) )
                if(m>0)
                    val = 0;
                else
                    val = 1;
                end
            else % elseif( (mod(n,4)==2) || (mod(n,4)==3) )
                if(m<0)
                    val = 0;
                else
                    val = 1;
                end
            end
            Znoll = n*(n+1)/2+abs(m)+val; %--index of Zernike Noll indices used
            
            indsZnoll(zz) = Znoll;
            zz = zz+1;
        end
    end
end


ZmapCubeInit = falco_gen_norm_zernike_maps(Nbeam,centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.


% for zz=1:Nmodes %floor(maxRadialOrder/2)+1
%     figure(20); imagesc(ZmapCubeInit(:,:,zz),1.5*[-1 1]); axis xy equal tight; colorbar;
%     title(sprintf('Mode %d',zz),'Fontsize',20); drawnow;
%     pause(0.1);
% end


end
%% Weed out modes if mirror x-axis of y-axis symmetry is required
if(strcmpi(symmaxis,'y')) %--Symmetry about y-axis 
    sumVals = zeros(Nmodes,1);
    for zz=1:Nmodes
        if(strcmpi(centering,'pixel'))
            Zmap = ZmapCubeInit(2:end,2:end,zz);
            Zmap = mask(2:end,2:end).*Zmap;
        else
            Zmap = ZmapCubeInit(:,:,zz);
            Zmap = mask.*Zmap;
        end
        sumVals(zz) = sum(sum(abs(Zmap+fliplr(Zmap))/2));
    end
    ZmapCube = ZmapCubeInit(:,:,sumVals>=1e-9);
    
elseif(strcmpi(symmaxis,'x')) %--Symmetry about x-axis 
    sumVals = zeros(Nmodes,1);
    for zz=1:Nmodes
        if(strcmpi(centering,'pixel'))
            Zmap = ZmapCubeInit(2:end,2:end,zz);
            Zmap = mask(2:end,2:end).*Zmap;
        else
            Zmap = ZmapCubeInit(:,:,zz);
            Zmap = mask.*Zmap;
        end
        sumVals(zz) = sum(sum(abs(Zmap+flipud(Zmap))/2));
    end
    ZmapCube = ZmapCubeInit(:,:,sumVals>=1e-9);
    
elseif(strcmpi(symmaxis,'none'))
    ZmapCube = ZmapCubeInit;
    
else
    error('Invalid input [%s] for symmetry in 3-fold Zernike generator function.',symmaxis)
end
    
Nz = size(ZmapCube,3);
fprintf('%d Zernike modes kept.\n',Nz);


% rmsVec = zeros(Nz,1);
% for zi=1:Nz
%     temp = ZmapCube(:,:,zi);
%     rmsVec(zi) = rms(temp(mask));
% end
% figure(32); plot(rmsVec);
% 
% for zi=1:size(ZmapCube,3) 
%     figure(40); imagesc(mask.*ZmapCube(:,:,zi),1.5*[-1 1]); axis xy equal tight; colorbar;
%     title(sprintf('Mode %d',zi),'Fontsize',20); drawnow;
%     pause(0.1);
% end


end %--END OF FUNCTION




