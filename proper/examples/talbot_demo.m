%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


% Demo script for the Talbot Effect example from the Proper Manual
%

% 2005 Feb     jek  created IDL routine
% 2017 Apr 13  gmg  Matlab translation
% 2017 Nov 13  gmg  Revised to match IDL variable names
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  diam   =    0.100d0;  % beam diameter (m)
  ftsz   =   16      ;  % size of font in plots
  ifig   =    0      ;  % figure number 
  n      =  128      ;  % number of pixels
  nseg   =    9      ;  % number of segments
  period =    0.040d0;  % period of cosine pattern (m)
  z      =    0.000d0;  % propagation distance (m)
  sa0    =  - 0.0013 ;  % scale of plot amplitude axis 0
  sa1    =    0.0013 ;  % scale of plot amplitude axis 1
  sp0    =  - 0.250d0;  % scale of plot phase axis 0
  sp1    =    0.250d0;  % scale of plot phase axis 1
  sx0    =    1      ;  % scale of plot X axis 0
  sx1    =   n       ;  % scale of plot X axis 1
  wavelength_microns = 0.500d0       ;  % wavelength (um)

  wavelength_m = wavelength_microns * 1.0d-6;

  talbot_length = 2.0d0 * period^2 / wavelength_m;      % Talbot length (m)
  delta_length  = talbot_length / (nseg - 1);

% Setup a 2 by "nseg" panel of plots
  figarray(nseg, 2);

  for i = 1 : nseg
    ov   = struct('diam', diam, 'dist', z, 'period', period);
    wavefront = prop_run('talbot', wavelength_microns, n, 'PASSVALUE', ov);

% Extract a horizontal cross-section of array
    wavefront = wavefront(fix(n / 2) + 1, :);

    amp  = abs(wavefront);
    amp  = amp - mean(amp);
    pha  = phase(wavefront);
    pha  = pha - mean(pha);

    ifig = ifig + 1;
    figplace(ifig);
    clf;                        % clear figure
    axes('FontSize', ftsz);
    plot(amp);
    axis([sx0, sx1, sa0, sa1]);
    set(get(gcf, 'CurrentAxes'), 'FontSize', ftsz);

    ifig = ifig + 1;
    figplace(ifig);
    clf;                        % clear figure
    axes('FontSize', ftsz);
    plot(pha);
    axis([sx0, sx1, sp0, sp1]);
    set(get(gcf, 'CurrentAxes'), 'FontSize', ftsz);

    z    = z + delta_length;
  end
