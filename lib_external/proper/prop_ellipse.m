%   Copyright 2016, 2017, 2020 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation of original version by Gary Gutt


function image = prop_ellipse(bm, rx, ry, varargin)
%        image = prop_ellipse(bm, rx, ry, varargin)
% Return an image containing an antialiased, filled ellipse
%
% Outputs:
%   image  = aperture mask containing antialiased filled ellipse
%
% Required Inputs:
% bm   = beam structure
% rx   = radius along x (meters unless norm = 1, then fraction of beam radius)
% ry   = radius along y (meters unless norm = 1, then fraction of beam radius)
%
% Optional inputs:
% 'xc'                = center of ellipse relative to wf center X
%                       (meters unless 'norm')
% 'yc'                = center of ellipse relative to wf center Y
%                       (meters unless 'norm')
% 'dark'              : if present, draw a dark ellipse (0 inside, 1 outside)
%                       (default is opposite way)
% 'norm'              = if present, radii and center coordinates are normalized
%                       to beam radius, otherwise they are in meters
% 'rotation'          = Degrees to rotate the ellipse about the center.
% 'subsample'         = Factor to subsample pixels in each dimension along ellipse 
%                       edge for antialiasing. Default is 11.  Must be odd-valued integer. 
%
% 2005 Feb     jek  created idl routine
% 2008 Sep     jek  Fixed bug that was causing a dark line when the
%                   ellipse edge was within 1/75 of the inner pixel boundary.
% 2014 May 13  gmg  Matlab translation
% 2015 Nov 17  gmg  Fixed bug in xmin, xmax to fill in ellipse correctly.
% 2016 Jun 14  gmg  Fixed bug for case where cx is outside of grid
% 2017 Mar 01  gmg  Revised for keyword/value for optional inputs
% 2020 Jan 27  jek  Added rotation (prototype algorithm suggested by A.J. Riggs (JPL)
%                   and subsampling keyword parameters; fixed 'norm' option
%                   so that simply specifying 'norm' will set it (cannot
%                   set it to zero in call)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    propcommon

    xc = 0.0;       % center of ellipse
    yc = 0.0;
    dark = 0;       % if set, have zero-filled interior
    norm = 0;       % radii, centers are normalized to beam radius 
    rotation = 0;   % rotation angle in degrees
    nsub = double(antialias_subsampling);      % subsampling factor used for antialiasing
    
    icav = 0;       % index in cell array varargin
    while icav < size(varargin, 2)
        icav = icav + 1;
        switch lower(varargin{icav})
            case {'cx', 'xc'}
                icav = icav + 1;
                xc   = varargin{icav};
            case {'cy', 'yc'}
                icav = icav + 1;
                yc   = varargin{icav};
            case {'dark'}
                dark = 1;
            case {'norm'}
                norm = 1;
            case {'rot','rotation'}
                icav = icav + 1;
                rotation = varargin{icav};
            otherwise
                error('prop_ellipse: Unknown keyword: %s\n', varargin{icav});
        end
    end
    
    dx = prop_get_sampling(bm);
    beamrad_pix = prop_get_beamradius(bm) / dx;
    [ny, nx] = size(bm.wf);
    n = nx;   
    xcenter_pix = floor(n/2) + 1;
    ycenter_pix = floor(n/2) + 1;
    
    if norm == 1
        xcenter_pix = xcenter_pix + xc * beamrad_pix;
        ycenter_pix = ycenter_pix + yc * beamrad_pix;
        xrad_pix = rx * beamrad_pix;                    % radius in X (pixels)
        yrad_pix = ry * beamrad_pix;                    % radius in Y (pixels)
    else
        xcenter_pix  = xcenter_pix + xc / dx;           % center X (pixels)
        ycenter_pix  = ycenter_pix + yc / dx;           % center Y (pixels)
        xrad_pix = rx / dx;                             % radius in X (pixels)
        yrad_pix = ry / dx;                             % radius in Y (pixels)
    end

    % rotate coordinates defining box containing unrotated ellipse

    sint = sind( rotation );
    cost = cosd( rotation );

    xp = [-xrad_pix xrad_pix xrad_pix -xrad_pix];
    yp = [-yrad_pix -yrad_pix yrad_pix yrad_pix];
    xbox = xp * cost - yp * sint + xcenter_pix;
    ybox = xp * sint + yp * cost + ycenter_pix;

    minx_pix = round(min(xbox)) - 1;
    maxx_pix = round(max(xbox)) + 1;
    minx_pix(minx_pix < 1) = 1;
    maxx_pix(maxx_pix > n) = n;
    nx = maxx_pix - minx_pix + 1;

    miny_pix = round(min(ybox)) - 1;
    maxy_pix = round(max(ybox)) + 1;
    miny_pix(miny_pix < 1) = 1;
    maxy_pix(maxy_pix > n) = n;
    ny = maxy_pix - miny_pix + 1;
    
    % create & rotate coordinate arrays 

    x = (0:nx-1) - xcenter_pix + minx_pix;
    x = repmat( x, ny, 1 );
    y = (0:ny-1) - ycenter_pix + miny_pix;
    y = repmat( reshape(y,ny,1), 1, nx );

    xr = (x * cost - y * sint) / xrad_pix;
    yr = (x * sint + y * cost) / yrad_pix;
    r = sqrt(xr.*xr + yr.*yr);
    drx = abs(r(1,2) - r(1,1));
    dry = abs(r(2,1) - r(1,1));

    delx = 1.0 / xrad_pix;
    dely = 1.0 / yrad_pix;
    drx = delx * cost - dely * sint;
    dry = delx * sint + dely * cost;

    dr = abs(drx);
    dr(dr < abs(dry)) = abs(dry);

    % find pixels along edge of ellipse
    
    mask = repmat(-1,ny,nx) .* (1 - (r > (1+dr)));
    mask(r <= (1-dr)) = 1;
    [wy, wx] = find(mask == -1);
    nw = length(wx);

    % subpixellate edge pixels and compute fractional coverage
    
    nsubpix = nsub * nsub;
    subpix_x_array = ((0:nsub-1) - floor(nsub/2)) / nsub + minx_pix - xcenter_pix;
    subpix_x_array = repmat( subpix_x_array, nsub, 1 );
    subpix_y_array = ((0:nsub-1) - floor(nsub/2)) / nsub + miny_pix - ycenter_pix;
    subpix_y_array = repmat( reshape(subpix_y_array,nsub,1), 1, nsub );
    limit = 1.0 + 1e-10;
    
    for k = 1:nw
        xs = subpix_x_array + wx(k) - 1;
        ys = subpix_y_array + wy(k) - 1;
        x = (xs * cost - ys * sint) / xrad_pix;
        y = (xs * sint + ys * cost) / yrad_pix;
        mask(wy(k),wx(k)) = sum(sum((x.*x+y.*y) <= limit)) / nsubpix;
    end
    
    if dark ~= 0
        image = ones(n,n);
        image(miny_pix:maxy_pix,minx_pix:maxx_pix) = 1 - mask;
    else
        image = zeros(n,n);
        image(miny_pix:maxy_pix,minx_pix:maxx_pix) = mask;
    end

end
