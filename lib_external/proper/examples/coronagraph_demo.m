%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


% Test script for run_coronagraph_dm.m

% 2005 Feb     jek  created idl routine
% 2017 Feb 08  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
% 2017 Oct 11  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Define plot parameters
  ftsz =     16;                % size of font in plots
  global ifig;                  % figure number
  ifig =      0;                % figure number

  n      =  512;                % grid size
  nps    =  256;                % number of pixels in psf sample
  nps2   = fix(nps / 2.0d0);
  icx    = fix(n   / 2.0d0) + 1;% index of center of array X
  icy    = fix(n   / 2.0d0) + 1;% index of center of array Y
  ix1    = icx - nps2;          % psf sample min index X
  ix2    = ix1 + nps - 1;       % psf sample max index X
  iy1    = icy - nps2;          % psf sample min index Y
  iy2    = iy1 + nps - 1;       % psf sample max index Y
  ca0    =    0.000d0;          % colorbar minimum
  ca1    =    0.030d0;          % colorbar maximum
  lambda =    0.550d0;          % wavelength (um)

  figarray(4, 4);               % set up positions of figures in array

  optval.use_dm = 0;                    % deformable mirror: no
  optval.use_errors = 0;                % prop_psd_errormap: no
  optval.occulter_type = '8TH_ORDER';   % prop_8th_order_mask: yes
  psf1 = prop_run('run_coronagraph_dm', lambda, n, 'prm', optval);
  psfs = psf1(iy1 : iy2, ix1 : ix2);

% Plot Point Spread Function intensity for case with no errors
  ifig = ifig + 1;
  figplace(ifig);
  clf;
  axes('FontSize', ftsz);
  imagesc(psfs.^0.25d0);
  axis equal;                   % x-axis units = y-axis units
  axis tight;                   % set axis limits to range of data
  axis xy;                      % set y-axis to increase from bottom
  hc = colorbar('vert');
  set(hc, 'FontSize', ftsz);
  caxis([ca0, ca1]);
  colormap(gray);
  set(get(gcf, 'CurrentAxes'), 'FontSize', ftsz);
  tit1 = sprintf('PSF: no errors');
  title(tit1, 'FontSize', ftsz);

  optval.use_dm = 0;                    % deformable mirror: no
  optval.use_errors = 1;                % prop_psd_errormap: yes
  optval.occulter_type = '8TH_ORDER';   % prop_8th_order_mask: yes
  psf2 = prop_run('run_coronagraph_dm', lambda, n, 'prm', optval);
  psfs = psf2(iy1 : iy2, ix1 : ix2);

% Plot Point Spread Function intensity for case with errors, no DM
  ifig = ifig + 1;
  figplace(ifig);
  clf;
  axes('FontSize', ftsz);
  imagesc(psfs.^0.25d0);
  axis equal;                   % x-axis units = y-axis units
  axis tight;                   % set axis limits to range of data
  axis xy;                      % set y-axis to increase from bottom
  hc = colorbar('vert');
  set(hc, 'FontSize', ftsz);
  caxis([ca0, ca1]);
  colormap(gray);
  set(get(gcf, 'CurrentAxes'), 'FontSize', ftsz);
  tit2 = sprintf('PSF: with errors');
  title(tit2, 'FontSize', ftsz);

  optval.use_dm = 1;                    % deformable mirror: yes
  optval.use_errors = 1;                % prop_psd_errormap: yes
  optval.occulter_type = '8TH_ORDER';   % prop_8th_order_mask: yes
  psf3 = prop_run('run_coronagraph_dm', lambda, n, 'prm', optval);
  psfs = psf3(iy1 : iy2, ix1 : ix2);

% Plot Point Spread Function intensity for case with errors, DM correction
  ifig = ifig + 1;
  figplace(ifig);
  clf;
  axes('FontSize', ftsz);
  imagesc(psfs.^0.25d0);
  axis equal;                   % x-axis units = y-axis units
  axis tight;                   % set axis limits to range of data
  axis xy;                      % set y-axis to increase from bottom
  hc = colorbar('vert');
  set(hc, 'FontSize', ftsz);
  caxis([ca0, ca1]);
  colormap(gray);
  set(get(gcf, 'CurrentAxes'), 'FontSize', ftsz);
  tit3 = sprintf('PSF: DM corrected');
  title(tit3, 'FontSize', ftsz);

  fprintf(1, 'Maximum speckle flux / stellar flux :\n');
  fprintf(1, '   No wavefront errors   = %16.7e\n', max(max(psf1)));
  fprintf(1, '   With wavefront errors = %16.7e\n', max(max(psf2)));
  fprintf(1, '   With DM correction    = %16.7e\n', max(max(psf3)));
