function B = falco_spatial_bin(I, binSize)
%spatialBin Strict block-sum binning for image cubes.
%
%   B = spatialBin(I, BIN)
%     - I: 2-D image [H x W] or 3-D image cube [H x W x N] (N slices).
%     - BIN: scalar (square) or [rows cols] block size.
%     - Sums over non-overlapping BIN-sized blocks in the spatial dims only.
%     - Throws an error if H or W is not divisible by BIN.
%
%   Notes
%     * No padding or cropping; strict divisibility required.
%     * The 3rd dimension (slices) is preserved and not mixed.
%     * Computation done in double, then cast back to input class.
%
    arguments
        I
        binSize {mustBeNumeric, mustBeNonempty}
    end

    % Normalize bin size
    if isscalar(binSize)
        by = binSize; bx = binSize;
    else
        by = binSize(1); bx = binSize(2);
    end
    if any([by bx] < 1) || any(mod([by bx],1) ~= 0)
        error('binSize must be positive integers (scalar or [rows cols]).');
    end

    % Accept 2-D or 3-D (cube) only
    sz = size(I);
    if numel(sz) == 2
        sz(3) = 1;                % treat 2-D as single-slice cube
    elseif numel(sz) > 3
        error('Only 2-D images or 3-D image cubes [H x W x N] are supported.');
    end
    [H, W, N] = deal(sz(1), sz(2), sz(3));

    % Strict divisibility check
    if mod(H, by) ~= 0 || mod(W, bx) ~= 0
        error(['Image size (%d x %d) must be divisible by binSize [%d %d]. ', ...
               'Consider padding externally or choosing a different binSize.'], ...
               H, W, by, bx);
    end

    % Number of blocks
    Hblocks = H / by;
    Wblocks = W / bx;

    % Work in double for safe summation
    A = double(I);

    % Reshape to [by, Hblocks, bx, Wblocks, N] then sum over by/bx
    A = reshape(A, by, Hblocks, bx, Wblocks, N);
    A = sum(A, 1);     % sum over block rows
    A = sum(A, 3);     % sum over block cols
    B = squeeze(A);    % -> [Hblocks x Wblocks x N]

    % Cast back like input
    B = cast(B, 'like', I);

    % If original was 2-D, return 2-D
    if N == 1
        B = B(:, :);
    end
end
