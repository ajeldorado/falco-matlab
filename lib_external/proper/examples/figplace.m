%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function hndl = figplace(ifig)
% function hndl = figplace(ifig)
% place a figure in the array of figures created by figarray
% Output:
% hndl = figure handle
% Input:
% ifig = figure number

% 2015 Jan 07  gmg  new routine

  global nfig pstn
  if ifig > nfig
    error('FIGPLACE: figure number %2d > number of figures %2d\n', ifig, nfig);
  else
    if ishghandle(ifig)
      hndl = figure(ifig);
    else
      hndl = figure('OuterPosition', pstn(ifig, :));
    end
  end
end                             % function figplace
