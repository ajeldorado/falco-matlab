%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function wavefront = coronagraph(wavefront, f_lens, occulter_type, diam)
%        wavefront = coronagraph(wavefront, f_lens, occulter_type, diam)
% A coronagraph example from the Proper Manual
%
% Outputs:
% wavefront     = beam structure (output)
%
% Required inputs:
% wavefront     = beam structure (input)
% f_lens        = focal length (m)
% occulter_type = occulter type
%               = '8TH_ORDER', 'GAUSSIAN', or 'SOLID'
% diam          = diameter (m)

% 2005 Feb     jek  created idl routine
% 2017 Feb 08  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
% 2017 Oct 11  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Define plot parameters
  ftsz       =   16;            % size of font in plots
  global ifig      ;            % index of figure

  wavefront  = prop_lens(wavefront, f_lens, 'coronagraph imaging lens');
  wavefront  = prop_propagate(wavefront, f_lens, 'snm', 'occulter');

% Occulter sizes are specified here in units of lambda / diameter.
% Convert lambda / diameter to radians, then to meters.

  lambda     = prop_get_wavelength(wavefront);      % wavelength (m)
  occrad     =    4.0d0;                % occulter radius in lambda / diam
  occrad_rad = occrad * lambda / diam;              % occulter radius (radians)
  dx_m       = prop_get_sampling(wavefront);        % pixel spacing (m)
  dx_rad     = prop_get_sampling_radians(wavefront);% pixel spacing (radians)
  occrad_m   = occrad_rad * dx_m / dx_rad;          % occulter radius (m)

  switch occulter_type
    case '8TH_ORDER'
      wavefront = prop_8th_order_mask(wavefront, occrad, ...
                    'tmin', 0.0d0, 'tmax', 1.0d0, 'circ');
    case 'GAUSSIAN'
      r          = prop_radius(wavefront);
      h          = sqrt(-0.5d0 * occrad_m^2 / log(1.0d0 - sqrt(0.5d0)));
      gauss_spot = 1.0d0 - exp(-0.5d0 * (r / h).^2);
      wavefront  = prop_multiply(wavefront, gauss_spot);
    case 'SOLID'
      wavefront  = prop_circular_obscuration(wavefront, occrad_m);
  end

% Plot (beam amplitude).^0.5
  ifig = ifig + 1;
  figure(ifig);
  clf;
  axes('FontSize', ftsz);
  imagesc((prop_get_amplitude(wavefront)).^0.5d0);
  axis equal;                   % x-axis units = y-axis units
  axis tight;                   % set axis limits to range of data
  axis xy;                      % set y-axis to increase from bottom
  hc = colorbar('vert');
  set(hc, 'FontSize', ftsz);
  colormap(gray);
  set(get(gcf, 'CurrentAxes'), 'FontSize', ftsz);
  tit1 = sprintf('After Occulter');
  title(tit1, 'FontSize', ftsz);

  wavefront  = prop_propagate(wavefront, f_lens, 'snm', 'pupil reimaging lens');
  wavefront  = prop_lens(wavefront, f_lens, 'pupil reimaging lens');

  wavefront  = prop_propagate(wavefront, f_lens * 2.0d0, 'snm', 'lyot stop');

% Plot (beam amplitude)^0.2
  ifig = ifig + 1;
  figure(ifig);
  clf;
  axes('FontSize', ftsz);
  imagesc((prop_get_amplitude(wavefront)).^0.2d0);
  axis equal;                   % x-axis units = y-axis units
  axis tight;                   % set axis limits to range of data
  axis xy;                      % set y-axis to increase from bottom
  hc = colorbar('vert');
  set(hc, 'FontSize', ftsz);
  colormap(gray);
  set(get(gcf, 'CurrentAxes'), 'FontSize', ftsz);
  tit2 = sprintf('Before Lyot Stop');
  title(tit2, 'FontSize', ftsz);

  switch occulter_type
    case '8TH_ORDER'
      wavefront  = prop_circular_aperture(wavefront, 0.50d0, 'norm');
    case 'GAUSSIAN'
      wavefront  = prop_circular_aperture(wavefront, 0.25d0, 'norm');
    case 'SOLID'
      wavefront  = prop_circular_aperture(wavefront, 0.84d0, 'norm');
  end

  wavefront  = prop_propagate(wavefront, f_lens, 'snm', 'reimaging lens');
  wavefront  = prop_lens(wavefront, f_lens, 'reimaging lens');

  wavefront  = prop_propagate(wavefront, f_lens, 'snm', 'final focus');

end                     % function coronagraph
