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
%   'sectors', and 'staircase', 'wrapped'.
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

    % OPTIONAL INPUTS
    centering = 'pixel';  %--Default to pixel centering
    xOffset = 0;
    yOffset = 0;
    clocking = 0; % [degrees]
    if(isfield(inputs,'centering')); centering = inputs.centering; end % 'pixel' or 'interpixel'
    if(isfield(inputs, 'xOffset')); xOffset = inputs.xOffset; end % [pixels]
    if(isfield(inputs, 'yOffset')); yOffset = inputs.yOffset; end % [pixels]
    if(isfield(inputs,'clocking')); clocking = inputs.clocking; end % [degrees]

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
            
            
            cent = [N/2+1, N/2+1];     % center for polar transform (center of vort in this case)
            rmax = min([size(vort)-max(cent), cent]);  % Maximum radius to which transform should be computed
            radpts = rmax;    % Number of points to interpolate on in radius
            angpts = 360;     % Number of points to interpolate on in angle
            [vortRad, rvec, qvec] = polarTransform(vort, cent, rmax, radpts, angpts, 'linear');  % Compute polar transform
%             figure(65); plot(qvec, mean(vortRad)/pi); xlim([0 2*pi]);LineWidth = 3;title('Vortex Phase Profile'); xlabel('\Theta / \pi'); ylabel('G_{0}/ \pi');


        case{'cos'}
            z_m = besselzero(0,1,1);
            z_m = z_m(end);
            vort = phaseScaleFac*z_m*cos(charge*THETA);
            mask = exp(1j*vort);
            
            cent = [N/2+1, N/2+1];     % center for polar transform (center of vort in this case)
            rmax = min([size(vort)-max(cent), cent]);  % Maximum radius to which transform should be computed
            radpts = rmax;    % Number of points to interpolate on in radius
            angpts = 360;     % Number of points to interpolate on in angle
            [vortRad, rvec, qvec] = polarTransform(vort, cent, rmax, radpts, angpts, 'linear');  % Compute polar transform
            figure(66); plot(qvec, mean(vortRad)/pi); xlim([0 2*pi]);LineWidth = 3;title('Cos Phase Profile'); xlabel('\Theta / \pi'); ylabel('G_{0}/ \pi');

             % Display vortex as a function of radius
             % Each row in the vortRad matrix represents a radius
             % Each col in the vortRad matrix represents an angle (coordinate angle in the pupil)
%            figure(); imagesc(rad2deg(qvec), rvec, vortRad);  colorbar; 


        case 'sectors'
            vort = phaseScaleFac*pi/2*sign(cos(charge*THETA));
            mask = exp(1j*vort);
            
            cent = [N/2+1, N/2+1];     % center for polar transform (center of vort in this case)
            rmax = min([size(vort)-max(cent), cent]);  % Maximum radius to which transform should be computed
            radpts = rmax;    % Number of points to interpolate on in radius
            angpts = 360;     % Number of points to interpolate on in angle
            [vortRad, rvec, qvec] = polarTransform(vort, cent, rmax, radpts, angpts, 'linear');  % Compute polar transform
            figure(67); plot(qvec, mean(vortRad)/pi); xlim([0 2*pi]);LineWidth = 3;title('Sectors Phase Profile'); xlabel('\Theta / \pi'); ylabel('G_{0}/ \pi');
                    
        case 'frenchwrapped'
            
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

            %  Display
            % Plot the vortex 
            figure(); imagesc(fancyVort); axis image; colorbar; title(' French Wrapped Vortex Phase Map');

            % Compute the radial profile
            cent = [N/2+1, N/2+1];     % center for polar transform (center of vort in this case)
            rmax = min([size(fancyVort)-max(cent), cent]);  % Maximum radius to which transform should be computed
            radpts = rmax;    % Number of points to interpolate on in radius
            angpts = 360;     % Number of points to interpolate on in angle
            [vortRad, rvec, qvec] = polarTransform(fancyVort, cent, rmax, radpts, angpts, 'linear');  % Compute polar transform

            % Since the vortex phase is radially symmetric, plot a radial average   
            figure(); plot(qvec, mean(vortRad)/pi); xlim([0 2*pi]);LineWidth = 3;title('French Wrapped Phase Profile'); xlabel('\Theta / \pi'); ylabel('G_{0}/ \pi');
            mask = exp(phaseScaleFac*1j*fancyVort);
            
            
        case 'classicalwrapped'
            
            vort = phaseScaleFac*charge*rem(THETA,pi./4);
            mask = exp(1j*vort);
            
                        
            cent = [N/2+1, N/2+1];     % center for polar transform (center of vort in this case)
            rmax = min([size(vort)-max(cent), cent]);  % Maximum radius to which transform should be computed
            radpts = rmax;    % Number of points to interpolate on in radius
            angpts = 360;     % Number of points to interpolate on in angle
            [vortRad, rvec, qvec] = polarTransform(vort, cent, rmax, radpts, angpts, 'linear');  % Compute polar transform
            figure(78); plot(qvec, mean(vortRad)/pi); xlim([0 2*pi]);LineWidth = 3;title('Classical Wrapped Phase Profile'); xlabel('\Theta / \pi'); ylabel('G_{0}/ \pi');

            
%             coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
%             %
%             fancyVort = 0.* coords.THETA;
%             
%             domain = (coords.THETA > 0) & (coords.THETA < pi);
%             vara = floor(coords.THETA(domain)./ (pi./8));
%             fancyVort(domain) = 8*coords.THETA(domain) - vara*2*pi;
%             domain = (coords.THETA > -pi) & (coords.THETA < 0);
%             vara = floor(coords.THETA(domain)./ (pi./8));
%             fancyVort(domain) = 8*coords.THETA(domain) - vara*2*pi;
%             domain = (coords.THETA > 0) & (coords.THETA < pi/8);
%             fancyVort(domain) = 8*coords.THETA(domain);
% 
%             domain = (coords.THETA > pi/2) & (coords.THETA < 5*pi/8);
%             fancyVort(domain) = 8*coords.THETA(domain) - 4*pi;
%             domain = (coords.THETA > 5*pi/8) & (coords.THETA < pi);
%             fancyVort(domain) = 8*coords.THETA(domain) - 6*pi;
%             domain = (coords.THETA > -pi) & (coords.THETA < -5*pi/8);
%             fancyVort(domain) = 8*(coords.THETA(domain)+pi);
%             domain = (coords.THETA > -pi+3*pi/8) & (coords.THETA < -pi/2);
%             fancyVort(domain) = 8*(coords.THETA(domain)+pi) - 2*pi;
%             domain = (coords.THETA > -pi/2) & (coords.THETA < -3*pi/8);
%             fancyVort(domain) = 8*(coords.THETA(domain)+pi) - 4*pi;
%             domain = (coords.THETA > -3*pi/8) & (coords.THETA < 0);
%             fancyVort(domain) = 8*(coords.THETA(domain)+pi) - 6*pi;

%             %  Display
%             % Plot the vortex 
%             figure(); imagesc(fancyVort); axis image; colorbar; title('Classical Wrapped Vortex Phase Map');
% 
%             % Compute the radial profile
%             cent = [N/2+1, N/2+1];     % center for polar transform (center of vort in this case)
%             rmax = min([size(fancyVort)-max(cent), cent]);  % Maximum radius to which transform should be computed
%             radpts = rmax;    % Number of points to interpolate on in radius
%             angpts = 360;     % Number of points to interpolate on in angle
%             [vortRad, rvec, qvec] = polarTransform(fancyVort, cent, rmax, radpts, angpts, 'linear');  % Compute polar transform
% 
%             % Since the vortex phase is radially symmetric, plot a radial average   
%             figure(); plot(qvec, mean(vortRad)/pi); xlim([0 2*pi]);LineWidth = 3;title('Classical Wrapped Phase Profile'); xlabel('\Theta / \pi'); ylabel('G_{0}/ \pi');
%             mask = exp(phaseScaleFac*1j*fancyVort);

        case 'staircase'
            Nsteps = inputs.Nsteps;
            vort = phaseScaleFac*ceil(mod((THETA+pi)/(2*pi)*charge, 1)*Nsteps)/Nsteps*2*pi;
            mask = exp(1j*vort);
%             figure(68); plot(THETA/pi,mask); axis xy equal tight; title('staircase phase function');drawnow;
            % mask = exp(1j*falco_gen_spiral_staircase(inputs)); % Uses antialiased edges.
            
            
            cent = [N/2+1, N/2+1];     % center for polar transform (center of vort in this case)
            rmax = min([size(vort)-max(cent), cent]);  % Maximum radius to which transform should be computed
            radpts = rmax;    % Number of points to interpolate on in radius
            angpts = 360;     % Number of points to interpolate on in angle
            [vortRad, rvec, qvec] = polarTransform(vort, cent, rmax, radpts, angpts, 'linear');  % Compute polar transform
            figure(69); plot(qvec, mean(vortRad)/pi); xlim([0 2*pi]);LineWidth = 3;title('Staircase Phase Profile'); xlabel('\Theta / \pi'); ylabel('G_{0}/ \pi');


        otherwise
            validOptions = "Valid options are 'vortex', 'cos', 'sectors', 'staircase' and 'wrapped'.";
            error('%s is not a valid option for inputs.type. \n%s', inputs.type, validOptions)
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


    function x=besselzero(n,k,kind)
        k3=3*k;
        x=zeros(k3,1);
        for j=1:k3
            % Initial guess of zeros 
            x0=1+sqrt(2)+(j-1)*pi+n+n^0.4;
            % Do Halley's method
            x(j)=findzero(n,x0,kind);
            if x(j)==inf
                error('Bad guess.');
            end
        end
        x=sort(x);
        dx=[1;abs(diff(x))];
        x=x(dx>1e-8);
        x=x(1:k);
        function x=findzero(n,x0,kind)
        n1=n+1;     n2=n*n;
        % Tolerance
        tol=1e-12;
        % Maximum number of times to iterate
        MAXIT=100;
        % Initial error
        err=1;
        iter=0;
        while abs(err)>tol & iter<MAXIT
            switch kind
                case 1
                    a=besselj(n,x0);    
                    b=besselj(n1,x0);   
                case 2
                    a=bessely(n,x0);
                    b=bessely(n1,x0);
            end
            x02=x0*x0;
            err=2*a*x0*(n*a-b*x0)/(2*b*b*x02-a*b*x0*(4*n+1)+(n*n1+x02)*a*a);
            x=x0-err;
            x0=x;
            iter=iter+1;
        end
        if iter>MAXIT-1
            warning('Failed to converge to within tolerance. ',...
                    'Try a different initial guess');
            x=inf;    
        end
        end
    end
        
end

