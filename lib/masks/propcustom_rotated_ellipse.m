

function apm = prop_rotated_ellipse(bm, rx, ry, varargin)
%        apm = prop_ellipse(bm, rx, ry, varargin)
% Return an image containing an antialiased, filled ellipse
%
% Outputs:
% apm  = aperture mask containing antialiased filled ellipse
%
% Required Inputs:
% bm   = beam structure
% rx   = radius along x (meters unless norm = 1, then fraction of beam radius)
% ry   = radius along y (meters unless norm = 1, then fraction of beam radius)
%
% Optional inputs:
% 'xc'                = center of ellipse relative to wf center X
%                       (m unless 'norm')
% 'yc'                = center of ellipse relative to wf center Y
%                       (m unless 'norm')
% 'dark'              : draw a dark circle (0 inside, 1 outside)
%                       (default is opposite way)
% 'norm'              = 1 radii and center coordinates are normalized
%                       to beam radius
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    % Set default values of input parameters
    cx   = 0.0;                           % center of rectangle X
    cy   = 0.0;                           % center of rectangle Y
    dark = 0;
    norm = 0;
    rot = 0;

    % Set values of internal parameters
    %   del  = 0.0000001;
    dx   = prop_get_sampling(bm);         % spacing between points in x (m)
    dy   = prop_get_sampling(bm);         % spacing between points in y (m)
    [ny, nx] = size(bm.wf);               % number of pixels in wavefront array
    prx  = prop_get_beamradius(bm) / dx;  % beam radius x in pixels
    pry  = prop_get_beamradius(bm) / dy;  % beam radius y in pixels

  icav = 0;                     % index in cell array varargin
    while icav < size(varargin, 2)
        icav = icav + 1;
        switch lower(varargin{icav})
            case {'cx', 'xc'}
                icav = icav + 1;
                cx   = varargin{icav};
            case {'cy', 'yc'}
                icav = icav + 1;
                cy   = varargin{icav};
            case {'dark'}
                dark = 1;
            case {'norm'}
                icav = icav + 1;
                norm = varargin{icav};
            case {'rot','rotation'}
                icav = icav + 1;
                rot   = varargin{icav};
            otherwise
                error('prop_rotated_ellipse: Unknown keyword: %s\n', varargin{icav});
        end
    end

    Nbeam = 2.0*prop_get_beamradius(bm) / dx; %inputs.Nbeam; % max aperture radius in samples
    Narray = nx; %inputs.Narray;% Number of samples across in square output array
    if norm == 1
        cpx  = floor(nx / 2 + 1) + cx * prx;
        cpy  = floor(ny / 2 + 1) + cy * pry;
        radx = rx * prx;                    % radius in X (pixels)
        rady = ry * pry;                    % radius in Y (pixels)
        radiusX = radx/Nbeam; %inputs.radiusX; % x-radius of ellipse [pupil diameters]
        radiusY = rady/Nbeam; %inputs.radiusY; % y-radius of ellipse [pupil diameters]
        xShear = cx * prx/Nbeam;
        yShear = cy * pry/Nbeam;
    else
        cpx  = floor(nx / 2 + 1) + cx / dx;   % center X (pixels)
        cpy  = floor(ny / 2 + 1) + cy / dy;   % center Y (pixels)
        radx = rx / dx;                       % radius in X (pixels)
        rady = ry / dy;                       % radius in Y (pixels)
        radiusX = radx/Nbeam; %inputs.radiusX; % x-radius of ellipse [pupil diameters]
        radiusY = rady/Nbeam; %inputs.radiusY; % y-radius of ellipse [pupil diameters]
        xShear = cx / (2*prop_get_beamradius(bm));
        yShear = cy / (2*prop_get_beamradius(bm));
    end


    clockingDegrees = rot; %inputs.clockingDegrees; % clocking of the pupil [degrees]
    

    magFac = 1.;
    
%     centering = 'pixel';
%     switch centering
%         case 'interpixel'
%             x = (-(Narray-1)/2:(Narray-1)/2)/Nbeam;
%         otherwise
%             x = (-Narray/2:(Narray/2-1))/Nbeam;
%     end
    
    x = (-Narray/2:(Narray/2-1))/Nbeam;
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

    halfWindowWidth = max(abs([RHO(2,1)-RHO(1,1), RHO(1,2) - RHO(1,1)]));%dx;
    apm = -1*ones(size(RHO));
    apm(abs(RHO) < radius - halfWindowWidth) = 1;
    apm(abs(RHO) > radius + halfWindowWidth) = 0;
    grayInds = find(apm==-1);
%     fprintf('Number of grayscale points = %d\n', length(grayInds));

    upsampleFactor = 101;%100;
    dxUp = dx/upsampleFactor;
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
        apm(grayInds(iInterior)) = pixelValue;
    %     figure(8); imagesc(xUp+xCenter, xUp+yCenter, subpixel); axis xy equal tight; colorbar; drawnow;

    end
    
    if dark == 1
        apm  = 1.0 - apm;
    end
    figure(6); imagesc(x, x, apm); axis xy equal tight; colorbar; drawnow;
%     keyboard
end % END OF FUNCTION