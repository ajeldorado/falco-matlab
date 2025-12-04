% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Store normalization values in the out structure.
%

function out = falco_store_normalization_data(mp, out, Itr)

    out.I00compactZeros = mp.I00compactZeros; % Star normalization for mp.dm1.V = mp.dm2.V = zeros
    out.I00fullZeros = mp.I00fullZeros; % Star normalization for mp.dm1.V = mp.dm2.V = zeros
    out.I00compactHist(:, Itr) = mp.Fend.compact.I00(:);
    out.I00fullHist(:, :, Itr) = mp.Fend.full.I00;
    out.I00compactRatioHist(:, Itr) = mp.Fend.compact.I00ratio;
    out.I00fullRatioHist(:, :, Itr) = mp.Fend.full.I00ratio;

end
