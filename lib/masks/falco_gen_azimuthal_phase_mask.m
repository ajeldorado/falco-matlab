% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% function OUT = propcustom_mft_PtoFtoP(IN, FPM, charge, apRad, inVal, outVal, useGPU ) 
%
% Generate a 
%
% FPM must be an array with dimensions MxM, where M = lambdaOverD*(beam diameter)
%
% INPUTS
% ------
% inputs: structure of inputs parameters
%   - inputs.type: type of mask. Valid options are 'vortex', 'cos',
%   'sectors','staircase','sawtooth','galicher8','wrapped6','dzpm','custom'
%   - inputs.charge: number of 2-pi phase progressions over the 360
%     degrees of the mask)
%   - inputs.N: width and height of the output array
%   - inputs.phaseScaleFac: Factor to apply uniformly to the phase.
%                           Used to add chromaticity.
%   - inputs.Nsteps: (required for 'staircase') number of steps per
%                    2*pi radians in the staircase
%   - inputs.clocking: (optional) clocking of the phase mask in degrees
%   - inputs.centering: (optional) array centering of the mask
%   - inputs.xOffset: (optional) x-offset of the mask center from the
%                     array center in units of pixels
%   - inputs.yOffset: (optional) y-offset of the mask center from the
%                     array center in units of pixels
%   - inputs.wav: (required for 'custom') index of fpm cube with lam depth
%
% OUTPUTS
% -------
% OUT: 2-D, NxN array with complex transmission of the mask 

function mask = falco_gen_azimuthal_phase_mask(inputs)

    % Required inputs
    maskType = inputs.type; 
    charge = inputs.charge; 
    N = inputs.N;
    phaseScaleFac = inputs.phaseScaleFac;
    
    
    %ADDITIONAL TEMPORARY PARAMETERS--needs to be properly integrated...
%     
%     P1 = 300; %mp.P1.full.Nbeam
%     pixPerLamD = 8;%N/P1;
    

    % OPTIONAL INPUTS
    centering = 'pixel';  %--Default to pixel centering
    xOffset = 0;
    yOffset = 0;
    clocking = 0; % [degrees]
    res = 0;
    if(isfield(inputs,'centering')); centering = inputs.centering; end % 'pixel' or 'interpixel'
    if(isfield(inputs, 'xOffset')); xOffset = inputs.xOffset; end % [pixels]
    if(isfield(inputs, 'yOffset')); yOffset = inputs.yOffset; end % [pixels]
    if(isfield(inputs,'clocking')); clocking = inputs.clocking; end % [degrees]
    if(isfield(inputs,'res'));res = inputs.res; end %[m]
    if(isfield(inputs,'roddierradius')); roddierradius = inputs.roddierradius; end % [lambda/D]
    if(isfield(inputs,'roddierphase')); roddierphase = inputs.roddierphase; end % [waves]
    % Input checks
    Check.scalar_integer(charge);
    
    % make coordinate system 
    if strcmpi(centering, 'pixel')
        xs = -N/2:(N/2-1);
    elseif strcmpi(centering, 'interpixel')
        xs = -(N-1)/2:(N-1)/2;
    else
        error('Invalid option for inputs.centering')
    end
    [X, Y] = meshgrid(xs);
    
    X = X - xOffset;
    Y = Y - yOffset;
    
    % rotate all points before computing THETA
    xyAll = [reshape(X, [1, N*N]); reshape(Y, [1, N*N])];
    rotMat = [cosd(clocking), sind(clocking); -sind(clocking), cosd(clocking)];
    xyAll = rotMat * xyAll;
    X = reshape(xyAll(1, :), [N, N]);
    Y = reshape(xyAll(2, :), [N, N]);

    [THETA, ~] = cart2pol(X, Y);

    % make mask 
    switch lower(maskType)

        case 'vortex'
            vort = phaseScaleFac*charge*THETA;
            mask = exp(1j*vort);
            
%             figure(); imagesc(vort); axis image; axis off; colorbar('Ticks',charge*pi.*[-0.99,-.5,0,0.5,1],...
%          'TickLabels',{num2str(-charge)+"\pi",num2str(-charge/2)+"\pi",'0',num2str(charge/2)+"\pi",num2str(charge)+"\pi"},'FontSize',20); title('Classical Vortex Phase Map','FontSize',20);

            
%           %FIGURE MAKING

%             set(groot,'defaulttextinterpreter','latex');
%             set(groot,'defaultLegendInterpreter','latex');
%             set(groot,'defaultAxesTickLabelInterpreter','latex');  
%             
%             figure(); imagesc(vort); axis image; axis off; set(gcf,'color','w'); set(gca,'fontsize',40); set(gcf,'position',[30,30,600,400])
%             h = colorbar('Ticks',charge*pi.*[-0.99,-.5,0,0.5,1],'TickLabels',{'$-6\pi$','$-3\pi$','$0$','$3\pi$','$6\pi$'},'TickLabelInterpreter', 'latex'); 
%             set(get(h,'label'),'string','\phi');
            %             title('Classical Vortex Phase Map','FontSize',20);

        case{'cos'}
            z_m = besselzero(0,1);
            z_m = z_m(end);
            vort = phaseScaleFac*z_m*cos(charge*THETA);
            mask = exp(1j*vort);

        case 'sectors'
            % Generate 1D phase profile 
            Nfac = 100;
            Nsamps1D = Nfac*N; 
            q = linspace(-pi,pi,Nsamps1D);
            fpmPhz1D = pi/2*sign(cos(charge*q));
            
            % Filter out unwanted modes  
            mask1D = exp(1j*phaseScaleFac*fpmPhz1D);
            mask1D_FT = fft(mask1D);
            filter = abs(mask1D_FT) > max(abs(mask1D_FT))/Nfac;
            mask1D_FT = mask1D_FT.*filter; 
            mask1D = ifft(mask1D_FT);

            % Convert to 2D mask 
            mask = interp1(q,mask1D,THETA,'linear');
%             mask = abs(mask).*exp(1j*phaseScaleFac*angle(mask));
            
            disp('generating sector mask')
            
        case 'staircase'
            Nsteps = inputs.Nsteps;

            % Generate 1D phase profile 
            Nfac = 1000;
            Nsamps1D = Nfac*N; 
            q = linspace(-pi,pi,Nsamps1D);
            fpmPhz1D = ceil(mod((q+pi)/(2*pi)*charge, 1)*Nsteps)/Nsteps*2*pi-pi;
            
            % Filter out unwanted modes  
            mask1D = exp(1j*phaseScaleFac*fpmPhz1D);
            mask1D_FT = fft(mask1D);
            filter = abs(mask1D_FT) > max(abs(mask1D_FT))/Nfac;
            mask1D_FT = mask1D_FT.*filter; 
            mask1D = ifft(mask1D_FT);

            % Convert to 2D mask 
            mask = interp1(q,mask1D,THETA,'linear');
%             figure;imagesc(angle(mask));colorbar; colormap(hsv);caxis([-pi pi]);set(gca,'ydir','normal')
%             mask = abs(mask).*exp(1j*phaseScaleFac*angle(mask));
            
            
            vort = phaseScaleFac*ceil(mod((THETA+pi)/(2*pi)*charge, 1)*Nsteps)/Nsteps*2*pi;
            
%             figure(); imagesc(vort); axis image; axis off; colorbar('Ticks',pi.*[0,0.5,1,1.5,1.99],...
%          'TickLabels',{0,"\pi/2","\pi","3\pi/2","2\pi"},'FontSize',20);title(' Staircase Vortex Phase Map','FontSize',20);


%           %FIGURE MAKING

%             set(groot,'defaulttextinterpreter','latex');
%             set(groot,'defaultLegendInterpreter','latex');
%             set(groot,'defaultAxesTickLabelInterpreter','latex');  
%             
%             figure(); imagesc(vort); axis image; axis off; set(gcf,'color','w');  set(gca,'fontsize',40); set(gcf,'position',[30,30,600,400])
%             h = colorbar('Ticks',pi.*[0,0.5,1,1.5,1.99],'TickLabels',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'},'TickLabelInterpreter', 'latex'); 
%             set(get(h,'label'),'string','\phi');
%             
        case 'sawtooth'
            
%             vort = phaseScaleFac*charge*rem(THETA,pi./4);
            coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
            vort = 0.* coords.THETA;
            domain = (coords.THETA >= 0);
            vort(domain) = charge*rem(coords.THETA(domain),2*pi./charge);
            domain = (coords.THETA >= -pi) & (coords.THETA < 0);
            vort(domain) = charge*rem((coords.THETA(domain)+pi),2*pi./charge);
            mask = exp(phaseScaleFac*1j*vort);

%             figure(); imagesc(vort); axis image; axis off; colorbar('Ticks',pi.*[0,0.5,1,1.5,1.99],...
%          'TickLabels',{0,"\pi/2","\pi","3\pi/2","2\pi"},'FontSize',20);title(' Sawtooth Vortex Phase Map','FontSize',20);

            %FIGURE MAKING
% 
%             set(groot,'defaulttextinterpreter','latex');
%             set(groot,'defaultLegendInterpreter','latex');
%             set(groot,'defaultAxesTickLabelInterpreter','latex');  
%             
%             figure(); imagesc(vort); axis image; axis off; set(gcf,'color','w');  set(gca,'fontsize',40); set(gcf,'position',[30,30,600,400])
%             h = colorbar('Ticks',pi.*[0,0.5,1,1.5,1.99],'TickLabels',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'},'TickLabelInterpreter', 'latex'); 
%             set(get(h,'label'),'string','\phi');


        case 'wrapped8'
            
            coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
            %
            fancyVort = 0.* coords.THETA;
            domain = (coords.THETA > 0) & (coords.THETA < 3*pi/8);
            fancyVort(domain) = 8*coords.THETA(domain);
            domain = (coords.THETA > 3*pi/8) & (coords.THETA < pi/2);
            fancyVort(domain) = 8*coords.THETA(domain) - 2*pi;
            domain = (coords.THETA > pi/2) & (coords.THETA < 5*pi/8);
            fancyVort(domain) = 8*coords.THETA(domain) - 4*pi;
            domain = (coords.THETA > 5*pi/8) & (coords.THETA < pi);
            fancyVort(domain) = 8*coords.THETA(domain) - 6*pi;
            domain = (coords.THETA > -pi) & (coords.THETA < -5*pi/8);
            fancyVort(domain) = 8*(coords.THETA(domain)+pi);
            domain = (coords.THETA > -pi+3*pi/8) & (coords.THETA < -pi/2);
            fancyVort(domain) = 8*(coords.THETA(domain)+pi) - 2*pi;
            domain = (coords.THETA > -pi/2) & (coords.THETA < -3*pi/8);
            fancyVort(domain) = 8*(coords.THETA(domain)+pi) - 4*pi;
            domain = (coords.THETA > -3*pi/8) & (coords.THETA < 0);
            fancyVort(domain) = 8*(coords.THETA(domain)+pi) - 6*pi;
            mask = exp(phaseScaleFac*1j*fancyVort);

%             figure(); imagesc(fancyVort); axis image; axis off; colorbar('Ticks',pi.*[-0.99,0,1,2,2.99],...
%          'TickLabels',{"-\pi",0,"\pi","2\pi","3\pi"},'FontSize',20);title(' Galicher Wrapped Vortex Phase Map','FontSize',20);

     
        case 'wrapped6'
            
            coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
            %
            fancyVort = 0.* coords.THETA;            
            domain = (coords.THETA > 0) & (coords.THETA < 1.84799568);
            fancyVort(domain) = 6*coords.THETA(domain);
            domain = (coords.THETA > 1.84799568) & (coords.THETA < 2.52559409);
            fancyVort(domain) = 6*coords.THETA(domain) - 2*pi;
            domain = (coords.THETA > 2.52559409) & (coords.THETA <= pi);
            fancyVort(domain) = 6*coords.THETA(domain) - 4*pi;
            domain = (coords.THETA > -pi+2.52559409) & (coords.THETA < 0);
            fancyVort(domain) = 6*(coords.THETA(domain)+pi)-4*pi;
            domain = (coords.THETA > -pi+1.84799568) & (coords.THETA < -pi+2.52559409);
            fancyVort(domain) = 6*(coords.THETA(domain)+pi)-2*pi;
            domain = (coords.THETA >= -pi) & (coords.THETA < -pi+1.84799568);
            fancyVort(domain) = 6*(coords.THETA(domain)+pi);
            mask = exp(phaseScaleFac*1j*fancyVort);
            
%             figure(); imagesc(fancyVort); axis image; axis off; colorbar('Ticks',pi.*[0,1,2,2.99],...
%          'TickLabels',{0,"\pi","2\pi","3\pi"},'FontSize',20);title('Wrapped Vortex Phase Map','FontSize',20);

        case 'dzpm'
            % dual zone phase mask
            if~res
               error('Error. For radial FPMs, the resolution must be specified.')
            end

            coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
            vort = 0.* coords.THETA;
            domain = (coords.THETA >= 0);
            vort(domain) = charge*rem(coords.THETA(domain),2*pi./charge);
            domain = (coords.THETA >= -pi) & (coords.THETA < 0);
            vort(domain) = charge*rem((coords.THETA(domain)+pi),2*pi./charge);
            
            R1 = (coords.RHO <= 0.515*res);
            vort(R1) =vort(R1) + 0.47*2*pi;
            R2 = (coords.RHO > 0.515*res) & (coords.RHO <= 0.705*res);
            vort(R2) =vort(R2) + 0.92*2*pi;
            
            mask = exp(phaseScaleFac*1j*vort);

            
        case 'roddier'
            %roddier + sawtooth
            if ~res
               error('Error. For radial FPMs, the resolution must be specified.')
            end

            if ~isfield(inputs, 'roddierradius')
                error("inputs.roddierradius must be defined for this mask case.")
            end

            if ~isfield(inputs, 'roddierphase')
                error("inputs.roddierphase must be defined for this mask case.")
            end

            coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
            vort = 0.* coords.THETA;
            domain = (coords.THETA >= 0);
            vort(domain) = charge*rem(coords.THETA(domain),2*pi./charge);
            domain = (coords.THETA >= -pi) & (coords.THETA < 0);
            vort(domain) = charge*rem((coords.THETA(domain)+pi),2*pi./charge);
            
            R1 = (coords.RHO <= roddierradius*res*phaseScaleFac);
            vort(R1) =vort(R1) + roddierphase*2*pi;
            
            mask = exp(phaseScaleFac*1j*vort);

%             figure(); imagesc(vort); axis image; colorbar('Ticks',pi.*[0,0.5,1,1.5,1.99],...
%          'TickLabels',{0,"\pi/2","\pi","3\pi/2","2\pi"},'FontSize',20);title('Roddier Phase Map','FontSize',20);

        case 'just_dimple'
            % a roddier dimple without any azimuthal structure
            if~res
                error('Error. For radial FPMs, the resolution must be specified.')
            end

            if ~isfield(inputs, 'roddierradius')
                error("inputs.roddierradius must be defined for this mask case.")
            end

            if ~isfield(inputs, 'roddierphase')
                error("inputs.roddierphase must be defined for this mask case.")
            end

            coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
            dimple = 0.* coords.THETA;
            domain = (coords.THETA >= 0);
            
            R1 = (coords.RHO <= roddierradius*res*phaseScaleFac);
            dimple(R1) =dimple(R1) + roddierphase*2*pi;
            
            mask = exp(phaseScaleFac*-1j*dimple);

%             figure(); imagesc(vort); axis image; colorbar('Ticks',pi.*[0,0.5,1,1.5,1.99],...
%          'TickLabels',{0,"\pi/2","\pi","3\pi/2","2\pi"},'FontSize',20);title('Roddier Phase Map','FontSize',20);


        case 'custom'
            % MONOCHROMATIC metasurface vortex
            
            %read in Metasurface vortcube, lams
%             load lorenzovortex6cube.mat
%             lams   = [3.4,3.6,3.8,4.,4.2];
            i_lam = 3; %find(round((lams(ceil(end/2)))./lams,4) == round(phaseScaleFac,4));
            
            
%             mask = vortcube(:,:,i_lam);
            
            design="sawtooth";
            trans  = "var"; %"const";
%             
%             % i_lam  = 2;%2;%3;%4;%5
% %             vortcube = [];%zeros(4);
%             
% 
% %             for i_lam = 1:length(lams)
% 
                Npad   = N; %4096; %4096; %8192; % %2^13;
                charge = 6;
                theta  = angle(falco_gen_vortex_mask(1,Npad)); % gen_v_mask  with charge 1 means having an array with the azi angle for each pxl
                z1 = 2.40483; % wiki Bessel zeros
                z2 = 5.52008; % using z2 or z3 results in more than 2pi phase coverage needed
                z3 = 8.65373; %
                if design=="sawtooth", [pss,Tss,sss] = b30_library_siz_auto(angle(exp(1i.*charge.*theta)),i_lam); end % VORTEX: add wavelength dependent phase. 
                if design=="cosine", [pss,Tss,sss] = b30_library_siz_auto(angle(exp(1i.*z1.*cos(charge.*theta)))+z1-pi,i_lam); end % COSINE: add wavelength dependent phase. 
                vortex1 = exp(1i.*pss); % phase only
                if trans=="const", vortex = vortex1; end
                if trans=="var", vortex = vortex1.*sqrt(Tss); end % add phase dependent transmission values
% 
% %                 vortcube = cat(3,vortcube,vortex); %(:,:,i_lam);
% 
% %             end
%             
            mask = vortex;
            
            
            
%             figure(); imagesc(angle(mask)); axis image; title('Custom Phase Map','FontSize',20);
%             disp('read metasurface map successfully');

            
        otherwise
            validOptions = "Valid options are 'vortex', 'cos', 'sectors', 'staircase', 'sawtooth', 'wrapped8', 'wrapped6', 'dzpm', 'roddier', and 'just_dimple'.";
            error('%s is not a valid option for inputs.type. \n%s', inputs.type, validOptions)
    end

    % Apply phase scaling factor to account for chromaticity 
%     mask = abs(mask).*exp(1j*phaseScaleFac*angle(mask));
    
end

function coords = generateCoordinates( N )
    %[ X,Y,THETA,RHO,xvals,yvals ] = generateCoordinates( N )
    %   Generates sample centered coordinate system (both cartesian and polar)

        % Create coordinate system 
        [X,Y] = meshgrid(-N/2:N/2-1);
        [THETA,RHO] = cart2pol(X,Y);
        xvals = X(1,:);yvals = Y(:,1);

        coords.N = N;
        coords.X = X;
        coords.Y = Y;
        coords.THETA = THETA;
        coords.RHO = RHO;
        coords.xvals = xvals;
        coords.yvals = yvals;
end

function [pt,rvec,qvec] = polarTransform(input_image, center_vec, rmax, numRadPts, numAngles,method)
    %[pt,rvec,qvec] = polarTransform(input_image, center_vec, rmax, numRadPts, numAngles,method)
    %   Computes the polar transform of input_image.
    %
    %   Inputs:
    %       input_image: 2D array
    %       center_vec: coordinates of origin in the image
    %       rmax: max radius to compute polar transform 
    %       numRadPts: number of points in the radial direction
    %       numAngles: number of polar angles
    %       method: interpolation method. See interp2 docs
    %   Outputs:
    %       pt: 2D array with polar transform
    %       rvec: vector of radial positions 
    %       qvec: vector of polar angles 


        rvec = linspace(0,rmax,numRadPts);% make array of radial points
        qvec = linspace(0,2*pi-2*pi/numAngles,numAngles);% make array of azimuthal points
        radialComp = repmat(rvec', 1, numAngles);% 2D array of radial points 
        angleComp = repmat(qvec, numRadPts, 1);% 2D array of azimuthal points 
        [xComp,yComp] = pol2cart(angleComp,radialComp);% Transform desired polar coords into cartesian coords

        % Make image coordinates 
        [rows,cols] = size(input_image);

        xvals = (1:rows) - center_vec(1);
        yvals = (1:cols) - center_vec(2);

        % compute polar transform
        pt = interp2(xvals,yvals,double(input_image),xComp,yComp,method);
end

% function x = besselzero(n,k,kind)
%         k3=3*k;
%         x=zeros(k3,1);
%         for j=1:k3
%             % Initial guess of zeros 
%             x0=1+sqrt(2)+(j-1)*pi+n+n^0.4;
%             % Do Halley's method
%             x(j)=findzero(n,x0,kind);
%             if x(j)==inf
%                 error('Bad guess.');
%             end
%         end
%         x=sort(x);
%         dx=[1;abs(diff(x))];
%         x=x(dx>1e-8);
%         x=x(1:k);
%         function x=findzero(n,x0,kind)
%         n1=n+1;     n2=n*n;
%         % Tolerance
%         tol=1e-12;
%         % Maximum number of times to iterate
%         MAXIT=100;
%         % Initial error
%         err=1;
%         iter=0;
%         while abs(err)>tol & iter<MAXIT
%             switch kind
%                 case 1
%                     a=besselj(n,x0);    
%                     b=besselj(n1,x0);   
%                 case 2
%                     a=bessely(n,x0);
%                     b=bessely(n1,x0);
%             end
%             x02=x0*x0;
%             err=2*a*x0*(n*a-b*x0)/(2*b*b*x02-a*b*x0*(4*n+1)+(n*n1+x02)*a*a);
%             x=x0-err;
%             x0=x;
%             iter=iter+1;
%         end
%         if iter>MAXIT-1
%             warning('Failed to converge to within tolerance. ',...
%                     'Try a different initial guess');
%             x=inf;    
%         end
%         end
% end
