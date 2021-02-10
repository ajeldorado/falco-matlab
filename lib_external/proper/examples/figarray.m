%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function pstn = figarray(nrow, ncol)
% function pstn = figarray(nrow, ncol)
% create position coordinates for an array of figures
% Output:
% pstn = array of coordinates of figures
% Inputs:
% nrow = number of rows
% ncol = number of columns

% 2015 Jan 07  gmg  new routine

  global nfig pstn
  gap  = 24;                            % gap on sides of each figure (pixels)
  mp   = get(0, 'MonitorPositions');
  nfig = ncol * nrow;                   % number of figures
  pstn = zeros(nfig, 4);                % figure positions (pixels)
  wdth = 560;                           % figure width (pixels)
  if ncol > 3
    spcx = (mp(1, 3) - wdth - gap) / (ncol - 1);% figure spacing x (pixels)
  else
    spcx = wdth + gap;                  % figure spacing x (pixels)
  end
  hght = 497;                           % figure height (pixels)
  if nrow > 2
    spcy = (mp(1, 4) - hght - gap) / (nrow - 1);% figure spacing y (pixels)
  else
    spcy = hght;                        % figure spacing y (pixels)
  end

  ifig = 0;                             % figure number
  for ir = 0 : nrow - 1
    for ic = 0 : ncol - 1
      ifig = ifig + 1;
      pstn(ifig, 1) = gap + ic * spcx;
      pstn(ifig, 2) = mp(1, 4) - hght - gap - ir * spcy;
      pstn(ifig, 3) = wdth;
      pstn(ifig, 4) = hght;
    end
  end
end                             % function figarray
