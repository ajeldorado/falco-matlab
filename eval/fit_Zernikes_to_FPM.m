
%--Library locations
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path


FPM = 1e9*fitsread('/Users/ajriggs/Repos/falco-matlab/data/brief/data_DM9_surf.fits');

FPMflat = fitsread('/Users/ajriggs/Repos/falco-matlab/data/brief/data_DM8_surf.fits');
FPMflat = FPMflat/max(FPMflat(:));
mask = FPMflat==1;

figure(8); imagesc(FPMflat); axis xy equal tight; colorbar; drawnow;
figure(9); imagesc(FPM); axis xy equal tight; colorbar; drawnow;
figure(10); imagesc(mask); axis xy equal tight; colorbar; drawnow;

%%
num_z = 250;

radius = 54;

[ z_coef, FPM_fit ] = prop_fit_zernikes( FPM, mask, radius, num_z);

figure(11); imagesc(FPM_fit.*mask); axis xy equal tight; colorbar; drawnow;
figure(12); imagesc((FPM_fit-FPM).*mask); axis xy equal tight; colorbar; drawnow;

figure(21); semilogy(1:num_z,abs(z_coef));

% [ z_coef, fitted_map ] =
% prop_fit_zernikes( input_map, mask, radius, num_z
% [, 'OBSCURATION_RATIO', value ] [, 'XC', xcenter ] [, 'YC', ycenter ] );