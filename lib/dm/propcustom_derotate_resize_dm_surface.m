%   Copyright 2016, 2017, 2019 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt
%   Modified from prop_dm by A.J. Riggs to fit the actuator commands to a given surface. 

function gridDerotAtActRes = propcustom_derotate_resize_dm_surface(surfaceToFit, dx, Nact, dmcx, dmcy, spcg, varargin)

% Derotate and resize a surface to the size and alignment of the actuator grid.
%
% Outputs:
% bm   = beam structure
% map  = returns DM surface (not wavefront) on same grid as bm (m)
%
% Required inputs:
% bm   = beam structure
% dmz0 = either:
%        2D array containing the surface piston of each DM actuator (m)
%        or the name of a 2D FITS (Flexible Image Transport System)
%        image file containing the above
% dmcx, dmcy = the location of the optical axis (center of the
%        wavefront) on the DM in actuator units (0 to nAct - 1).
%        The center of the first actuator is (0.0, 0.0).
% spcg = defines the spacing between actuators (m);
%        must not be used when n_act_acroos_pupil is specified.
%
% Optional inputs:
% 'xtilt'             = specify the rotation of the DM surface with
% 'ytilt'               respect to the wavefront plane, with the origin
% 'ztilt'               at the center of the wavefront.  The DM surface
%                       is interpolated and orthogonally projected onto
%                       the wavefront grid.  The coordinate system
%                       assumes that the wavefront and initial DM
%                       surface are in the X,Y plane with a lower left
%                       origin with Z towards the observer.  The
%                       rotations are left handed.  The default rotation
%                       order is X, Y, then Z unless the zyx switch is
%                       set. (default: 0, 0, 0)
% 'zyx'               : specifies the rotation order if two or more of
%                       xtilt, ytilt, or ztilt are specified.
%                       The default is X, Y, then Z rotations.
%                       Set for for Z, Y, then X rotations.
% 'inf_fn','inf_file' : specify a new influence function as a FITS file
%                       with the same header keywords as PROPER's default
%                       influence function. Needs these values in
%                       info.PrimaryData.Keywords:
%                       'P2PDX_M' % pixel width x (m)
%                       'P2PDY_M' % pixel width y (m)
%                       'C2CDX_M' % actuator pitch x (m)
%                       'C2CDY_M' % actuator pitch y (m)
% 'inf_sign'          : specifies the sign (+/-) of the influence function.
%                       Given as an option because the default influence
%                       function file is positive, but positive DM actuator
%                       commands make a negative deformation for Xinetics
%                       and BMC DMs.
%
% 2005 Feb     jek  created idl routine 
% 2007 Jun     jek  fixed bug in interpolation of smoothed DM surface &
%                   changed default from smoothed to non-smoothed
% 2014 Mar     jek  implemented faster method of creating DM surface
% 2014 Jun 26  gmg  Matlab translation
% 2015 Oct     jek  Implemented DM surface tilting, rotation
% 2015 Dec 04  gmg  Matlab translation
% 2017 Feb 15  gmg  Revised for keyword/value for optional inputs
% 2019 Feb 15  a r  Added two new keyword/value pairs as optional inputs:
%                   -Accept any influence function from a FITS file
%                   -Allow the sign of the influence function to be + or -
% 2019 Nov 20  a r  Changed to de-rotate and resize a surface to the size 
%                   of the DM command array. To get the voltages, call 
%                   falco_fit_dm using the result of this function.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dmz = eye(Nact); %--Only used for its size. The actual values are thrown away.
%   if ischar(dmz0)               % then open 2D FITS image file
%     dmz  = fitsread(dmz0);
%   else
%     dmz  = dmz0;
%   end

%   Fit  = 0;             % values in dmz0 are commanded actuator heights
%   nAct = 0;             % number of actuators across pupil
%   NoAp = 0;             % the DM pattern is added to the wavefront
  tlt  = zeros(1, 3);   % set default tilts to 0 degrees
  zyx  = 0;             % rotation order is X, Y, then Z rotations
  DMInfFuncFileName = 'influence_dm5v2.fits'; % default influence function, Xinetics
  sign_factor = 1;      % positive or negative sign multiplied with influence function

  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
%       case {'fit'}
%         Fit  = 1;       % values in dmz0 are desired surface heights
%       case {'noap', 'no_apply'}
%         NoAp = 1;       % the DM pattern is not added to the wavefront
%       case {'nact', 'n_act_across_pupil'}
%         icav = icav + 1;
%         nAct = varargin{icav};  % number of actuators across pupil
      case {'tltx', 'xtilt'}
        icav = icav + 1;
        tlt(1) = varargin{icav};% X rotation of the DM surface (deg)
      case {'tlty', 'ytilt'}
        icav = icav + 1;
        tlt(2) = varargin{icav};% Y rotation of the DM surface (deg)
      case {'tltz', 'ztilt'}
        icav = icav + 1;
        tlt(3) = varargin{icav};% Z rotation of the DM surface (deg)
      case {'xyz'}
        zyx  = 0;       % rotation order is X, Y, then Z rotations
      case {'zyx'}
        zyx  = 1;       % rotation order is Z, Y, then X rotations
      case{'inf_fn','inf_file'} % name of FITS file for influence function
        icav = icav + 1;
        DMInfFuncFileName = varargin{icav};
      case{'inf_sign'} % + or - sign for influence function
        icav = icav + 1;
        sign_char = varargin{icav};
        switch lower(sign_char(1))
            case{'+','p'} % positive or plus sign
                sign_factor = 1;
            case{'-','n','m'} % negative or minus sign
                sign_factor = -1;
            otherwise
                error('Chosen inf_sign value not allowed.')
        end
            
      otherwise
        error('Unknown keyword: %s\n', varargin{icav});
    end
  end

    % Read the influence function data from the specified FITS file
    info = fitsinfo(DMInfFuncFileName);
    inf  = fitsread(DMInfFuncFileName);
    
    %--Multiply by +1 or -1
    inf = sign_factor*inf; 

    [ldef, idef] = ismember('NAXIS1' , info.PrimaryData.Keywords(:, 1));
    ifnx = info.PrimaryData.Keywords{idef, 2};    % inf func number of pixels x

    [ldef, idef] = ismember('NAXIS2' , info.PrimaryData.Keywords(:, 1));
    ifny = info.PrimaryData.Keywords{idef, 2};    % inf func number of pixels y

    [ldef, idef] = ismember('P2PDX_M', info.PrimaryData.Keywords(:, 1));
    ifdx = info.PrimaryData.Keywords{idef, 2};    % inf func spacing x (m)

    [ldef, idef] = ismember('P2PDY_M', info.PrimaryData.Keywords(:, 1));
    ifdy = info.PrimaryData.Keywords{idef, 2};    % inf func spacing y (m)

    [ldef, idef] = ismember('C2CDX_M', info.PrimaryData.Keywords(:, 1));
    acdx = info.PrimaryData.Keywords{idef, 2};    % actuator spacing x (m)

    [ldef, idef] = ismember('C2CDY_M', info.PrimaryData.Keywords(:, 1));
    acdy = info.PrimaryData.Keywords{idef, 2};    % actuator spacing y (m)

  ifcx = floor(ifnx / 2) + 1;           % inf func center pixel x
  ifcy = floor(ifny / 2) + 1;           % inf func center pixel y
  raix = round(acdx / ifdx);            % ratio actuator / inf func dx
  raiy = round(acdy / ifdy);            % ratio actuator / inf func dy

%   if spcg ~=  0.0 & nAct ~=  0
%     error('Proper:PROP_DM', ...
%           'User cannot specify both actuator spacing and number of actuators across pupil.\n');
%   end
% 
%   if (spcg ==  0.0 & nAct ==  0)
%     error('Proper:PROP_DM', ...
%           'User must specify either actuator spacing or number of actuators across pupil.\n');
%   end

% Set the real DM actuator spacing
%   if nAct ~= 0
%     actd = 2.0 * prop_get_beamradius(bm) / nAct;
%   else
    actd = spcg;
%   end

% Scale the influence function sampling to the specified DM actuator spacing
  ifdx = ifdx * actd / acdx;            % actual inf func spacing x (m)
  ifdy = ifdy * actd / acdy;            % actual inf func spacing y (m)

%   if Fit == 1
% % Then calculate the actuator positions to fit the given DM surface
% % Create the meshgrids for the influence function inf
%     [cfx, cfy] = meshgrid([ 1 : ifnx], [ 1 : ifny]);
% % Create the meshgrids for the interpolated influence function infk
%     [ckx, cky] = meshgrid([-2 :  2  ], [-2 :  2  ]);
%     ckx  = ckx * actd;
%     cky  = cky * actd;
% % Calculate the interpolated influence function infk
%     infk = interp2(cfx, cfy, inf, ckx/ifdx + ifcx, cky/ifdy + ifcy, 'cubic');
% % Calculate the actuator positions to fit the given DM surface
%     dmcz = prop_fit_dm(dmz, infk);
%   else
% % Use the given actuator positions
%     dmcz = dmz;
%   end
  
  % Use the given actuator positions
    dmcz = dmz;

  [acny, acnx] = size(dmz);             % number of actuators in y and x
% Create subsampled DM grid
  mrgx = raix * 9;                      % margin on each end of each row
  mrgy = raiy * 9;                      % margin on each end of each column
  dmnx = raix * acnx + mrgx * 2;        % number of DM grid points x
  dmny = raiy * acny + mrgy * 2;        % number of DM grid points y
  lgx1 = fix(raix/2) + mrgx + 1;        % index of  1st actuator center x
  lgy1 = fix(raiy/2) + mrgy + 1;        % index of  1st actuator center y
  lgx2 = lgx1 + raix * (acnx-1);        % index of last actuator center x
  lgy2 = lgy1 + raiy * (acny-1);        % index of last actuator center y
  dmg0  = zeros(dmny, dmnx);             % initialize subsampled DM grid
% Fill the subsampled DM grid with the actuator amplitudes
  dmg0(lgy1 : raiy : lgy2, lgx1 : raix : lgx2) = dmcz;
% Do the convolution of the actuator amplitudes with the influence function
  dmg  = conv2(dmg0, inf, 'same');

  [ny, nx] = size(surfaceToFit);       % number of points in wavefront array
%   [ny, nx] = size(bm.wf);       % number of points in wavefront array

% 3D rotate DM grid and project orthogonally onto wavefront

% Calculate grid dimensions (pix) projected onto wavefront
  xdim = round(sqrt(2.0) * dmnx * ifdx / dx);
  if xdim > nx
    xdim = nx;
  end
  xd2  = fix(xdim / 2) + 1;
  ydim = round(sqrt(2.0) * dmny * ifdy / dx);
  if ydim > ny
    ydim = ny;
  end
  yd2  = fix(ydim / 2) + 1;

  cx   = ([1 : xdim] - xd2) * dx;
  cy   = ([1 : ydim] - yd2) * dx;
  [cxm, cym] = meshgrid(cx, cy);

  sa   = sind(tlt(1));
  ca   = cosd(tlt(1));
  sb   = sind(tlt(2));
  cb   = cosd(tlt(2));
  sg   = sind(tlt(3));
  cg   = cosd(tlt(3));

  if zyx == 0
    rotm = [               cb * cg,               -cb * sg,       sb, 0.0; ...
            ca * sg + sa * sb * cg, ca * cg - sa * sb * sg, -sa * cb, 0.0; ...
            sa * sg - ca * sb * cg, sa * cg + ca * sb * sg,  ca * cb, 0.0; ...
                               0.0,                    0.0,      0.0, 1.0];
  else
    rotm = [ cb * cg, sa * sb * cg - ca * sg, ca * sb * cg + sa * sg, 0.0; ...
             cb * sg, sa * sb * sg + ca * cg, ca * sb * sg - sa * cg, 0.0; ...
            -sb,      sa * cb,                ca * cb,                0.0; ...
                 0.0,                    0.0,                    0.0, 1.0];
  end
  
  %% Compute xdm0 and ydm0 for use in de-rotating the DM surface
  % Forward project a square
  edge = [-1.0, -1.0,  0.0,  0.0; ...
           1.0, -1.0,  0.0,  0.0; ...
           1.0,  1.0,  0.0,  0.0; ...
          -1.0,  1.0,  0.0,  0.0];
  xyzn = edge;% * rotm;   % had to reverse matrix mult. order to match IDL

% Determine backwards projection for screen-raster-to-DM-surface computation
% Had to reverse and increment indices to match IDL
  dxdx = (xyzn(1, 1) - xyzn(2, 1)) / (edge(1, 1) - edge(2, 1));
  dxdy = (xyzn(2, 1) - xyzn(3, 1)) / (edge(2, 2) - edge(3, 2));
  dydx = (xyzn(1, 2) - xyzn(2, 2)) / (edge(1, 1) - edge(2, 1));
  dydy = (xyzn(2, 2) - xyzn(3, 2)) / (edge(2, 2) - edge(3, 2));

  xs   = (cxm/dxdx - cym*dxdy / (dxdx*dydy)) / (1.0 - dydx*dxdy / (dxdx*dydy));
  ys   = (cym/dydy - cxm*dydx / (dxdx*dydy)) / (1.0 - dydx*dxdy / (dxdx*dydy));

  xdm0  = (xs + dmcx * actd) / ifdx + lgx1;
  ydm0  = (ys + dmcy * actd) / ifdy + lgy1;
  
  %% Compute xdm and ydm for use in de-rotating the DM surface
 
% Forward project a square
  edge = [-1.0, -1.0,  0.0,  0.0; ...
           1.0, -1.0,  0.0,  0.0; ...
           1.0,  1.0,  0.0,  0.0; ...
          -1.0,  1.0,  0.0,  0.0];
  xyzn = edge * rotm;   % had to reverse matrix mult. order to match IDL

% Determine backwards projection for screen-raster-to-DM-surface computation
% Had to reverse and increment indices to match IDL
  dxdx = (xyzn(1, 1) - xyzn(2, 1)) / (edge(1, 1) - edge(2, 1));
  dxdy = (xyzn(2, 1) - xyzn(3, 1)) / (edge(2, 2) - edge(3, 2));
  dydx = (xyzn(1, 2) - xyzn(2, 2)) / (edge(1, 1) - edge(2, 1));
  dydy = (xyzn(2, 2) - xyzn(3, 2)) / (edge(2, 2) - edge(3, 2));

  xs   = (cxm/dxdx - cym*dxdy / (dxdx*dydy)) / (1.0 - dydx*dxdy / (dxdx*dydy));
  ys   = (cym/dydy - cxm*dydx / (dxdx*dydy)) / (1.0 - dydx*dxdy / (dxdx*dydy));

  xdm  = (xs + dmcx * actd) / ifdx + lgx1;
  ydm  = (ys + dmcy * actd) / ifdy + lgy1;

  %% Derotate the DM surface
  grid = padOrCropEven(surfaceToFit,size(xdm,1));
  gridDerot = griddata(xdm, ydm, grid, xdm0, ydm0, 'cubic');
  gridDerot(isnan(gridDerot)) = 0;

%   figure(10); imagesc(gridDerot); axis xy equal tight; colorbar;
%   figure(11); imagesc(dmg(2:end,2:end)-rot90(dmg(2:end,2:end),2)); axis xy equal tight; colorbar;
%   figure(12); imagesc(gridDerot(2:end,2:end)-rot90(gridDerot(2:end,2:end),2)); axis xy equal tight; colorbar;

%% Resize and decimate the DM surface to get it at the same size as the DM actuator command array.
%  The result will be fed to falco_fit_dm_surf() for deconvolution with the
%  influence function.

  xOffsetInAct = ((Nact/2 - 1/2) - dmcx);
  yOffsetInAct = ((Nact/2 - 1/2) - dmcy);

  multipleOfCommandGrid = ceil_odd(spcg/dx);
  N1 = Nact*multipleOfCommandGrid;
  N2 = size(grid,1);
  xs1 = (-(N1-1)/2:(N1-1)/2)/N1 ; %--interpixel centered
  if(mod(N2,2)==0)
    xs2 = (-N2/2:(N2/2)-1)/N2*(N2*dx/(Nact*spcg));
  else
    xs2 = (-(N2-1)/2:(N2-1)/2)/N2*(N2*dx/(Nact*spcg));
  end
  [XS1,YS1] = meshgrid(xs1);
  [XS2,YS2] = meshgrid(xs2);

  gridDerotResize = interp2(XS2,YS2,gridDerot,XS1+xOffsetInAct/Nact,YS1+yOffsetInAct/Nact,'cubic', 0.0);
  
  gridDerotAtActRes = gridDerotResize(ceil(multipleOfCommandGrid/2):multipleOfCommandGrid:end,ceil(multipleOfCommandGrid/2):multipleOfCommandGrid:end); % decimate

%   dm.dm_spacing = spcg;
%   dm.dx_inf0 = 1e-4;%dx;%1e-4;
%   dm.inf0 = inf;
%   dm.Nact = Nact;
%   Vout = falco_fit_dm_surf(dm,gridDerotAtActRes);
%   
%   figure(19); imagesc(gridDerotResize); axis xy equal tight; colorbar;
%   figure(20); imagesc(dmz); axis xy equal tight; colorbar;
%   figure(21); imagesc(gridDerotAtActRes); axis xy equal tight; colorbar;
%   figure(22); imagesc(Vout); axis xy equal tight; colorbar;
% %   figure(23); imagesc(dmz0 - Vout); axis xy equal tight; colorbar;

end % END OF FUNCTION