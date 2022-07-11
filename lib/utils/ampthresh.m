% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Threshold a measured pupil to give a boolean mask of illuminated pixels.

function boolMask = ampthresh(pupilMap, varargin)

    %--Optional inputs
    if size(varargin, 2) == 1
        nBins = varargin{1};
    else
        nBins = 21;
    end

    pupilMap = abs(pupilMap); 

    % Input checks
    if ~isa(nBins, 'numeric')
        error('ampthresh:InputMustBeNumeric', 'Input must be numeric.')
    end
    if min(pupilMap(:)) == max(pupilMap(:))
        error('ampthresh:InputCannotBeUniform', 'Cannot threshold an array of uniform values.')
    end

    % use histogram of intensities to choose threshold
    [Icount, IbinEdges] = histcounts(pupilMap, nBins);

    % find the minima in the histogram
    bV = (Icount(2:end-1) <= Icount(1:end-2)) & (Icount(2:end-1) < Icount(3:end));

    % vector of intensity values at the center of each bin
    binCenter = 0.5*(IbinEdges(1:end-1) + IbinEdges(2:end));

    % List of bin values where the histogram has a miminum
    temp = binCenter(2:end-1);
    binVal = temp(bV);

    % Choose the first minimum as the threshold. this will be the lowest
    % with greater intensity must be "signal". If binVal is empty, the
    % histogram has no o minima except at an endpoint. This is a
    % failure, so return empty.
    if length(binVal) > 0
        thresh = binVal(1);
    else
        error('Failed to find a histogram minimum')
    end

    % Apply theshold to create boolean mask
    boolMask = pupilMap > thresh;

end

