% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% Function to create a rotated elliptical aperture.
%
% This function is needed because PROPER can generate only unrotated
% ellipses. The grayscale edges of the ellipse are generated by creating a
% binary-valued, higher-resolution sub-array of each pixel along the edge 
% of the ellipse.
%
%--Required Inputs in the structure "inputs"
% inputs.Nbeam % max aperture radius in samples
% inputs.Narray % Number of samples across in square output array
% inputs.radiusX % x-radius of ellipse [pupil diameters]
% inputs.radiusY % y-radius of ellipse [pupil diameters]
% inputs.clockingDegrees % clocking of the pupil [degrees]
%
%--Optional Inputs
% inputs.centering
% inputs.xShear % x-shear of ellipse [pupil diameters]
% inputs.yShear % y-shear of ellipse [pupil diameters]
% inputs.magFac % magnification factor


function pupil = falco_gen_ellipse(inputs)
    
    %--Required Inputs
    Nbeam = inputs.Nbeam; % max aperture radius in samples
    Narray = inputs.Narray;% Number of samples across in square output array
    radiusX = inputs.radiusX; % x-radius of ellipse [pupil diameters]
    radiusY = inputs.radiusY; % y-radius of ellipse [pupil diameters]
    clockingDegrees = inputs.clockingDegrees; % clocking of the pupil [degrees]
    
    %--Optional inputs
    centering = 'pixel';
    xShear = 0.;
    yShear = 0.;
    magFac = 1.;
    upsampleFactor = 100;
    if(isfield(inputs,'centering')); centering = inputs.centering; end
    if(isfield(inputs,'xShear')); xShear = inputs.xShear; end % x-shear of ellipse [pupil diameters]
    if(isfield(inputs,'yShear')); yShear = inputs.yShear; end % y-shear of ellipse [pupil diameters]
    if(isfield(inputs,'magFac')); magFac = inputs.magFac; end % magnification factor
    if(isfield(inputs,'upsampleFactor')); upsampleFactor = inputs.upsampleFactor; end % upsampling factor
    
    switch centering
        case 'interpixel'
            x = (-(Narray-1)/2:(Narray-1)/2)/Nbeam;
        otherwise
            x = (-Narray/2:(Narray/2-1))/Nbeam;
    end
    y = x;
    x = x - xShear;
    y = y - yShear;
    [X,Y] = meshgrid(x,y);
    dx = x(2) - x(1);
    radius = 0.5;

    RHO = 1/magFac*0.5*sqrt(...
        1/(radiusX)^2*(cosd(clockingDegrees)*X + sind(clockingDegrees)*Y).^2 + ...
        1/(radiusY)^2*(sind(clockingDegrees)*X - cosd(clockingDegrees)*Y).^2 ...
        );

    halfWindowWidth = max([radiusX/radiusY, radiusY/radiusX])*max(abs([RHO(2,1)-RHO(1,1), RHO(1,2) - RHO(1,1)]));
    pupil = -1*ones(size(RHO));
    pupil(abs(RHO) < radius - halfWindowWidth) = 1;
    pupil(abs(RHO) > radius + halfWindowWidth) = 0;
    grayInds = find(pupil==-1);
    % fprintf('Number of grayscale points = %d\n', length(grayInds));
    
%     rhoInterior = (pupil==-1);
%     xTransition = X(rhoInterior); 
%     yTransition = Y(rhoInterior); 

    dxUp = dx/upsampleFactor;
    % if mod(upsampleFactor, 2) == 1  % odd
    %     nhalf = floor((upsampleFactor-1)/2);
    %     xUp = (-nhalf:nhalf)*dxUp;
    % elseif mod(upsampleFactor, 2) == 0  % even
    %     nhalf = (upsampleFactor-1)/2;
    %     xUp = (-nhalf:nhalf)*dxUp;
    % end
    
    xUp = (-(upsampleFactor-1)/2:(upsampleFactor-1)/2)*dxUp;
    [Xup, Yup] = meshgrid(xUp);

    subpixel = zeros(upsampleFactor,upsampleFactor);

    for iInterior = 1:length(grayInds)

        subpixel = 0*subpixel;

        xCenter = X(grayInds(iInterior));
        yCenter = Y(grayInds(iInterior));
        RHOup = 0.5*sqrt(...
        1/(radiusX)^2*(cosd(clockingDegrees)*(Xup+xCenter) + sind(clockingDegrees)*(Yup+yCenter)).^2 + ...
        1/(radiusY)^2*(sind(clockingDegrees)*(Xup+xCenter) - cosd(clockingDegrees)*(Yup+yCenter)).^2 ...
        );

        subpixel((RHOup) <= radius) = 1;
        pixelValue = sum(subpixel(:))/upsampleFactor^2;
        pupil(grayInds(iInterior)) = pixelValue;
    %     figure(1); imagesc(xUp+xCenter, xUp+yCenter, subpixel); axis xy equal tight; colorbar; drawnow;

    end
    % figure(2); imagesc(x, x, pupil); axis xy equal tight; colorbar; drawnow;
    
end % END OF FUNCTION