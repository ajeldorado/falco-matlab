%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [bmo, mask] = prop_8th_order_mask(bmi, hwhm, varargin)
%        [bmo, mask] = prop_8th_order_mask(bmi, hwhm, varargin)
% Multiply the current wavefront by an 8th-order occulter.
% An 8th-order transmission mask is generated.  The mask can be linear
% (the default), circular, ar elliptical in shape.  The mask width is
% specified as the radius at which the transmission is 0.5 (along the
% X-axis for linear (unless yax is specified) and elliptical apertures).
%
% Outputs:
% bmo  = wavefront array structure, out
% mask = amplitude (not intensity) mask image array output
%
% Required inputs:
% bmi  = wavefront array structure, in
% hwhm = radius at which the mask transmission is 0.5
%        (along the direction of the X-axis for elliptical apertures).
%        By default, this is in units of lambda / D.
%        If 'METERS' is specified, this is in meters.
%
% Optional inputs:
% 'meters'            : if set, hwhm is in meters
% 'min_transmission'  = mask minimum intensity transmission; default is 0.0.
% 'max_transmission'  = mask maximum intensity transmission; default is 1.0.
% 'circular'          : if set, creates a circular occulter; default is linear.
% 'elliptical'        = creates an elliptical occulter with axis ratio
%                     = y_width / x_width.
% 'y_axis'            : if set, specifies that the linear occulter is to
%                       be drawn with transmission varying along the
%                       Y-axis rather than the X-axis (the default).
%                       Only valid for the linear occulter.

% 2005 Feb     jek  created idl routine
% 2016 Aug 08  gmg  Matlab translation
% 2017 Feb 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  circ =     0       ;  % linear occulter
  elra =     0.0d0   ;  % elliptical axis ratio
  tmax =     1.0d0   ;  % mask maximum intensity transmission
  tmin =     0.0d0   ;  % mask minimum intensity transmission
  yax  =     0       ;  % transmission varying along the X-axis

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'mtrs', 'meters'}
        hwhm = hwhm / (prop_get_fratio(bmi) * prop_get_wavelength(bmi));
      case {'tmin', 'min_transmission'}
        icav = icav + 1;
        tmin = varargin{icav};
      case {'tmax', 'max_transmission'}
        icav = icav + 1;
        tmax = varargin{icav};
      case {'circ', 'circular'}
        circ =     1 ;        % circular occulter
      case {'elra', 'elliptical'}
        icav = icav + 1;
        elra = varargin{icav};
      case {'yax', 'y_axis'}
        yax  =     1 ;        % transmission varying along the Y-axis
      otherwise
        error('prop_8th_order_mask: Unknown keyword: %s\n', varargin{icav});
    end
  end

  lnr  = (circ == 0) & (elra == 0);     % = 1 if linear

  dx   = prop_get_sampling(bmi);
  dy   = dx;
  [ny, nx] = size(bmi.wf);
  pare = 1.788 / hwhm;  % parameter epsilon in mask transmission eqn
  parl = 3d0;           % parameter l in mask transmission eqn
  parm = 1d0;           % parameter m in mask transmission eqn
  plml = (parl - parm) / parl;

% Compute beam pixel spacing in lambda / D units
  dxld = dx / (prop_get_fratio(bmi) * prop_get_wavelength(bmi));
  dyld = dy / (prop_get_fratio(bmi) * prop_get_wavelength(bmi));

  if lnr                % linear mask
    if yax              % along y-axis
% Create vectors of point coordinates y (L / D)
      vpy  = ([1 : ny] - (fix(ny / 2) + 1)) * dyld;
      mask = prop_sinc(vpy * pi * pare / parm).^parm * parm / parl ...
           - prop_sinc(vpy * pi * pare / parl).^parl + plml;
    else                % along x-axis
% Create vectors of point coordinates x (L / D)
      vpx  = ([1 : nx] - (fix(nx / 2) + 1)) * dxld;
      mask = prop_sinc(vpx * pi * pare / parm).^parm * parm / parl ...
           - prop_sinc(vpx * pi * pare / parl).^parl + plml;
    end

  elseif circ           % circular mask
% Create vectors of point coordinates (L / D)
    vpx  = ([1 : nx] - (fix(nx / 2) + 1)) * dxld;
    vpy  = ([1 : ny] - (fix(ny / 2) + 1)) * dyld;
    [apx, apy] = meshgrid(vpx, vpy);    % arrays of x, y (L / D)
    apr  = sqrt(apx.^2 + apy.^2);       % array of radii (L / D)
    mask = prop_sinc(apr * pi * pare / parm).^parm * parm / parl ...
         - prop_sinc(apr * pi * pare / parl).^parl + plml;

  else                  % elliptical mask
% Create vectors of point coordinates (L / D)
    vpx  = ([1 : nx] - (fix(nx / 2) + 1)) * dxld;
    vpy  = ([1 : ny] - (fix(ny / 2) + 1)) * dyld / elra;
    [apx, apy] = meshgrid(vpx, vpy);    % arrays of x, y (L / D)
    apr  = sqrt(apx.^2 + apy.^2);       % array of radii (L / D)
    mask = prop_sinc(apr * pi * pare / parm).^parm * parm / parl ...
         - prop_sinc(apr * pi * pare / parl).^parl + plml;
  end

% Renormalize mask
  mask = mask.^2;       % intensity
  mask = mask - min(min(mask));
  mask = mask / max(max(mask));
  mask = mask * (tmax - tmin) + tmin;
  mask = sqrt(mask);    % back to amplitude

% If linear, convert 1D mask into 2D mask
  if lnr
    if yax              % along y-axis
      [mskx, mask] = meshgrid(1 : nx, mask);
    else                % along x-axis
      [mask, msky] = meshgrid(mask, 1 : ny);
    end
  end

% Create output wavefront structure and apply amplitude mask
  bmo  = bmi;
  bmo.wf = bmi.wf .* prop_shift_center(mask);
end                     % function prop_8th_order_mask
