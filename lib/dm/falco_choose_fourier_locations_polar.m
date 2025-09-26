% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% [xiVec, etaVec] = falco_choose_fourier_locations_polar(freqMax, spacing, gridType, radiusInner, radiusOuter, angleOpen, clocking)
%
% Compute the (x,y) coordinates of all equally-spaced spatial frequencies in a given region.
%
% INPUTS
% freqMax: Maximum spatial frequency (along +/- x or y) to include.
% spacing: distance between spatial frequencies. Units of lambda/D.
% gridType: Type of grid. 'square' or 'hex'
% radiusInner: Inner radius to use. Units of lambda/D.
% radiusOuter: Outer radius to use. Units of lambda/D.
% angleOpen: Opening angle to use. Units of degrees.
% clocking: CCW rotation angle to use. Units of degrees.
% xiMin: (optional) Minimum allowed x-axis value (before clocking).
%   Use this to make a D-shape. Units of lambda/D.
%
% OUTPUT
% xiVec: vector of horizontal axis locations of the spatial frequencies
% etaVec: vector of vertical axis locations of the spatial frequencies
%
% NOTES
% Avoid zero frequency because that is just a piston on the DM, which
% should not do anything.


function [xiVec, etaVec] = falco_choose_fourier_locations_polar(freqMax, spacing, gridType, radiusInner, radiusOuter, angleOpen, clocking, varargin)

    %--Optional inputs
    if size(varargin, 2) == 1
        xiMin = varargin{1};
        Check.real_nonnegative_scalar(xiMin);
    else
        xiMin = [];
    end

    switch lower(gridType)

        case 'square'

            % 1st, Make grid of points within a giant square
            xisHalf = spacing/2:spacing:freqMax;
            xis = [-fliplr(xisHalf), xisHalf];
            [XIS, ETAS] = meshgrid(xis);


        case {'hex', 'hexagon', 'hexagonal'} % % TODO

            % 1st, Make vector of points within a giant hexagon
            Nrings = ceil(freqMax/(sqrt(3)/2)/spacing);
            XIS = [];
            ETAS = [];
            % row number (rowNum) is 1 for the center row and 2 is above it, etc.
            % Nacross is the total number of points across that row
            for rowNum = 1:Nrings
                Nacross = 2*Nrings - rowNum; % Number of actuators across at that row (for hex tiling in a hex shape)
                xiOffset = Nrings - (rowNum+1)/2; % x offset from origin

                xis = (0:Nacross-1).' - xiOffset; % xi values are 1 apart
                etas = sqrt(3)/2*(rowNum-1)*ones(Nacross,1); % same y-value for the entire row

                if rowNum == 1
                    XIS = [XIS; xis];
                    ETAS = [ETAS;etas];
                else
                    XIS = [XIS; xis; xis];
                    ETAS = [ETAS; etas; -etas]; % rows +/-n have +/- y coordinates
                end
            end
            % Rescale
            XIS = spacing * XIS;
            ETAS = spacing * ETAS;

        otherwise
            error("gridType must be 'square' or 'hex'.")

    end

    % 2nd, remove all points not in the region of interest
    [THETAS, RHOS] = cart2pol(XIS, ETAS);
%     THETAS_ROT = angle(exp(1j*(THETAS-deg2rad(clocking))));
%     THETAS_ROT = THETAS-deg2rad(clocking);
    radialMask = (RHOS >= radiusInner) & (RHOS <= radiusOuter);
    angleMask = abs(THETAS) <= deg2rad(angleOpen)/2;
    if ~isempty(xiMin)
        knifeEdgeMask = RHOS.*cos(THETAS) >= xiMin;
    else
        knifeEdgeMask = logical(ones(size(radialMask)));
    end
    softwareMask = radialMask & angleMask & knifeEdgeMask;

    xiVec = XIS(softwareMask);
    etaVec = ETAS(softwareMask);

    [THETAS, RHOS] = cart2pol(xiVec, etaVec);
    xiVec  = RHOS.*cos(THETAS - deg2rad(clocking));
    etaVec = RHOS.*sin(THETAS - deg2rad(clocking));
    

end %--END OF FUNCTION