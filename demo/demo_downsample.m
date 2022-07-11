% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------


%% Pixel centered array

clear

arrayIn = ones(3);
arrayIn = pad_crop(arrayIn, 9);
arrayIn = pad_crop(arrayIn, 15, 'extrapval', 1);
arrayIn = pad_crop(arrayIn, 21);

diam0 = 9;
diam1 = 3;
centering = 'pixel';

arrayOut = falco_filtered_downsample(arrayIn, diam1/diam0, centering);

arrayExpected = pad_crop(1, 3);
arrayExpected = pad_crop(arrayExpected, 5, 'extrapval', 1);
arrayExpected = pad_crop(arrayExpected, 7);

whos array*

figure(101); imagesc(arrayIn); axis xy equal tight; colorbar;
figure(102); imagesc(arrayOut); axis xy equal tight; colorbar;
figure(103); imagesc(arrayIn - fliplr(arrayIn)); axis xy equal tight; colorbar;
figure(104); imagesc(arrayOut - fliplr(arrayOut)); axis xy equal tight; colorbar;
figure(105); imagesc(arrayOut - arrayExpected); axis xy equal tight; colorbar;

%% Interpixel centered array

clear

arrayIn = ones(8);
arrayIn = pad_crop(arrayIn, 16);
arrayIn = pad_crop(arrayIn, 24, 'extrapval', 1);
arrayIn = pad_crop(arrayIn, 32);

diam0 = 8;
diam1 = 2;
centering = 'interpixel';

arrayOut = falco_filtered_downsample(arrayIn, diam1/diam0, centering);

whos array*

figure(101); imagesc(arrayIn); axis xy equal tight; colorbar;
figure(102); imagesc(arrayOut); axis xy equal tight; colorbar;
figure(103); imagesc(arrayIn - fliplr(arrayIn)); axis xy equal tight; colorbar;
figure(104); imagesc(arrayOut - fliplr(arrayOut)); axis xy equal tight; colorbar;
