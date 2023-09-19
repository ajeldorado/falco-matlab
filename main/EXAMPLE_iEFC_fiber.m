% Copyright 2018-2021 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to run wavefront control on a vortex coronagraph with a segmented
% input pupil.

clear

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
mp.flagSaveWS = true;  %--Whether to save out entire (large) workspace at the end of trial. Default is false


%% Step 2: Load default model parameters

EXAMPLE_defaults_iEFC_fiber
mp.estimator = 'iefc';
% mp.est.probe.axis = 'multi';

% % Fiber parameters
mp.flagFiber = true;
mp.flagLenslet = false;

mp.Fend.x_fiber = [6];%[5.3405 -2.6702 -2.6702]; %Fiber core center positions in lambda_0/D
mp.Fend.y_fiber = [0];%[0 4.625 -4.625];
mp.Fend.Nfiber = numel(mp.Fend.x_fiber);

mp.fiber.a = 0.507;%0.875;%0.66; %Radius of the fiber core in lambda_0/D
mp.fiber.a_phys = 1.75e-6; %Physical radius of the fiber core in meters
mp.fiber.NA = 0.12; %Numerical aperture of the fiber

%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%--Use just 1 wavelength for initial debugging/testing of code
mp.fracBW = 0.1;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 3;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

mp.Nitr = 5; %--Number of wavefront control iterations

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];
mp.runLabel_short = ['xSMF',num2str(mp.Fend.x_fiber),'_ySMF',num2str(mp.Fend.y_fiber),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100)];
outdir = strcat('/Users/jllopsay/Documents/MATLAB/hcstr/output/efcsmf_TTsens_2023/',mp.runLabel_short,'/');
mkdir(outdir)

%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);
mp.path.outdir = outdir;
%% Define the probes to use 

dmVsin = sum(mp.dm1.basisCube(:, :, 1:mp.dm1.NbasisModes/2), 3);
figure(31); imagesc(dmVsin); axis xy equal tight; colorbar; colormap parula;
dmVcos = sum(mp.dm1.basisCube(:, :, mp.dm1.NbasisModes/2+1:end), 3);
figure(32); imagesc(dmVcos); axis xy equal tight; colorbar; colormap parula;

mp.iefc.probeCube = zeros(mp.dm1.Nact, mp.dm1.Nact, 2);
mp.iefc.probeCube(:, :, 1) = dmVsin;
mp.iefc.probeCube(:, :, 2) = dmVcos;

mp.iefc.modeCoef = 1e-3; %--Gain coefficient to apply to the normalized DM basis sets for the empirical SCC calibration.
mp.iefc.probeCoef = 1e-3; %--Gain coefficient to apply to the stored probe commands used for IEFC state estimation.
mp.iefc.probeDM = 1; %--Which DM to use when probing for IEFC.
%% Verify visually that the Fourier modes fully cover the dark hole

freqMax = max([max(mp.dm1.fourier_basis_xis), max(mp.dm1.fourier_basis_etas)]);
figure(111);
imagesc(mp.Fend.xisDL, mp.Fend.etasDL, mp.Fend.corr.maskBool); colormap gray;
set(gca, 'Fontsize', 20);
set(gcf, 'Color', 'w');
hold on;
h111 = plot(mp.dm1.fourier_basis_xis , mp.dm1.fourier_basis_etas, 'or');
set(h111, 'MarkerFaceColor', 'r', 'MarkerSize', 5);
title('Overlay of Fourier Modes on the Dark Hole')
axis xy equal tight;
hold off;
drawnow;

%%
[mp, out] = falco_wfsc_loop(mp, out);

