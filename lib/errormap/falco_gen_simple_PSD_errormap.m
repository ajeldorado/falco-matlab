% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate the phase error of a mirror as 1/f^(alpha) noise.
%
% INPUTS
% ------
% N : the number of samples in each dimension. width of the output array
% alpha : the 1/f noise power
% mirror_figure : the peak-to-valley noise variation [arbitrary units]
%
% OUTPUTS
% -------
% out : the mirror phase noise in same units as mirror_figure
%
% Provided by Ruslan Belikov of NASA Ames Research Center.

function out = falco_gen_simple_PSD_errormap(N, alpha, mirror_figure)

    aberr = rand(N);
    FT = fft2(aberr);

    x = ((0:1:N-1)-N/2)' * ones(1, N);
    y = ones(N, 1) * ((0:1:N-1)-N/2);

    envelope = fftshift(1./(x.^2+y.^2).^(alpha/4));
    envelope(1, 1) = 0; % DC component undefined, set to 0

    FT = FT .* envelope;
    out = real(fftshift(ifft2(FT)));
    %out = out/(max(max(out)) - min(min(out)))*mirror_figure;
    out = out/sqrt(sum(sum(out.^2)))*mirror_figure*N;
    
end