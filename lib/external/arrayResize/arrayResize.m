function DataOut = arrayResize(DataIn,bin)
% Function to downsample an image stack based on averaging through a symmetric box filter.
% Usage: DataOut = arrayResize(DataIn,bin)
% DataIn = x*y matrix to be downsampled. Additional dimensions won't be affected.
% bin = size of the box filter. Downsampled image will consist of averages
% of squares of the size bin*bin. If DataIn can not be divided by bin, the
% downsampled image will not contain the lower and right edges of DataIn 
% that are above the highest divider.

dSize = size(DataIn); %size of input matrix
a = ['end-' num2str(rem(dSize(1),bin)-1) ':end'];
b = ['end-' num2str(rem(dSize(2),bin)-1) ':end'];

if rem(dSize(1),bin) > 0
    eval(['DataIn(' a repmat(',:',1,length(dSize)-1) ') = [];']); %delete overlapping edges
end
if rem(dSize(2),bin) > 0
    eval(['DataIn(:,' b repmat(',:',1,length(dSize)-2) ') = [];']); %delete overlapping edges
end

dSize = size(DataIn); %new size of input matrix
DataOut = mean(reshape(DataIn,[bin,prod(dSize(1:2))/bin,dSize(3:end)]),1,'native');
DataOut = reshape(DataOut,[dSize(1)/bin,dSize(2),dSize(3:end)]);

ind = 1:length(dSize);
DataOut = permute(DataOut,[ind(2) ind(1) ind(3:end)]); %permute matrix to average over 2nd dimension
dSize = size(DataOut); %size of permuted input matrix 
DataOut = mean(reshape(DataOut,[bin,prod(dSize(1:2))/bin,dSize(3:end)]),1,'native');
DataOut = reshape(DataOut,[dSize(1)/bin,dSize(2),dSize(3:end)]);
DataOut = permute(DataOut,[ind(2) ind(1) ind(3:end)]); %permute matrix back to original