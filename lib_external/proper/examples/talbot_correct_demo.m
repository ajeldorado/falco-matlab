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
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  diam =    0.100d0  ;  % beam diameter (m)
  ftsz =   16        ;  % size of font in plots
  ifig =    0        ;  % figure number 
  np   =  128        ;  % number of pixels
  nseg =    9        ;  % number of segments
  peri =    0.040d0  ;  % period of cosine pattern (m)
  pz   =    0.000d0  ;  % propagation distance (m)
  sa0  =  - 0.0013   ;  % scale of plot amplitude axis 0
  sa1  =    0.0013   ;  % scale of plot amplitude axis 1
  sp0  =  - 0.250d0  ;  % scale of plot phase axis 0
  sp1  =    0.250d0  ;  % scale of plot phase axis 1
  sx0  =    1        ;  % scale of plot X axis 0
  sx1  =   np        ;  % scale of plot X axis 1
  wlu  =    0.500d0  ;  % wavelength (um)

  wlm  = wlu * 1.0d-6;

  tal  = 2.0d0 * peri^2 / wlm;       % Talbot length (m)
  del  = tal / (nseg - 1);

  figarray(5, 4);               % set up positions of figures in array

  for it = 1 : nseg
    ov   = struct('diam', diam, 'dist', pz, 'peri', peri);
    wfa  = prop_run('talbot_correct', wlu, np, 'PASSVALUE', ov);

% Extract a horizontal cross-section of array
    wfax = wfa(fix(np / 2) + 1, :);

    amp  = abs(wfax);
    amp  = amp - mean(amp);
    pha  = phase(wfax);
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

    pz   = pz + del;
  end
