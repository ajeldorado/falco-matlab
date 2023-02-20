% Copyright 2018 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to run WFSC with a vortex using implicit EFC (IEFC) and a Fourier basis set.

% clear

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
addpath(genpath('/home/hcst/falco-matlab'));% savepath;
addpath('/home/hcst/HCST_scripts/EFC');
%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false

clearvars -except bench 

%% Step 2: Load default model parameters
lambda0 = 785*1e-9;

mp.flagFiber = false;
EXAMPLE_defaults_HCST_VVC % Re-use SCC example's configuration
mp.flagSim = false;      %--Simulation or not
mp.testbed = 'HCST';
mp.dm1.basisType = 'fourier';
mp.dm1.Nactbeam = 27;
mp.Fend.res = 8.7;

mp.Fend.sides = 'top'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'
mp.Fend.shape = 'D'; % 'D', 'circle'


mp.Nitr = 50;
%% HCST preparation
% Flags; all flags false --> regular EFC run
flag_svc = false;
efc_imageSharpening = false;
flagFieldStop = false;
mp.flag_lc = false;

% path2flatMap = [bench.info.HCST_DATA_DIR,'is_results/2022Oct13/'];
% path2flatMap =[bench.info.HCST_DATA_DIR,'pr_results/2022Nov23/'];
path2flatMap =[bench.info.HCST_DATA_DIR,'pr_results/2023Jan16/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'is_results/2022Oct13/'];
path2startingMap =[bench.info.HCST_DATA_DIR,'efc_results/2023Jan27/'];
% path2startingMap =[bench.info.HCST_DATA_DIR,'pr_results/2023Jan16/'];

NDfilter_FWpos = 2;
FWpos = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsbp = 1; % Number of sub-bandpasses, AKA number of wavelengths in full bandpass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FWpos~=NDfilter_FWpos
    NDfilter_cal = 28.6082 * 4.9367 * 4.4753 * 4.6956;%23.16; %>1
%     NDfilter_cal = 28.6082 * 4.9367 * 4.4753;% * 4.6956;%23.16; %>1
else
    NDfilter_cal = 1;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tint_offAxis = 1e-2*ones(Nsbp,1); %; % 1e-3; %Integration time for off-axis PSF
mp.tint_efc = 5e-0; %30; %1e-1; %30; %35e-0;%1e-1;%[2,1,1,0.5,1]*0.05; % for the control and evalution
mp.tint_est = 1e-2;%5e-1; % time of integration for the estimation/sensing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bench.info.sbp_texp = mp.tint_efc;
mp.tint = mp.tint_efc;
mp.peakPSFtint = tint_offAxis;
mp.NDfilter_cal = NDfilter_cal;

% Source parameters
bench.info.source = 'nkt';% 'nkt';%

if Nsbp>1
    fracBW = 0.10;  % bandwidth of experiment
    mp.sbp_centers = lambda0*linspace( 1-fracBW/2,1+fracBW/2,Nsbp);
    bench.info.sbp_width = 14e-9*ones(Nsbp,1);%[14,14,14,14,14]*1e-9;% 3e-9;%-  %--Width of each sub-bandpass on testbed (meters)
%     bench.info.sbp_width([4:6]) = 6e-9; %
else
    fracBW = 0.01; % bandwidth of experiment; when Nsbp is 1 --> narrowband or monochromatic  
    mp.sbp_centers = lambda0;
    bench.info.sbp_width = [14]*1e-9; % bandpass of the sub-bandpass [nm]
end

frameSize = 180;

datelabel = datestr(now,'yyyymmmddTHHMM');
label_progress = [num2str(Nsbp),'lams',num2str(round(1e9*lambda0)),'nm_BW',num2str(fracBW*100),'_VVC_',datelabel];

falcoPreparev2

mp.bench = bench;
sbp_width = 16e-9;
% bench.info.sbp_width = sbp_width ; %--Width of each sub-bandpass on testbed (meters)
bench.info.sbp_texp = mp.tint_efc;% Exposure time for each sub-bandpass (seconds)
bench.info.PSFpeaks = PSFpeaks;% counts per second 

%% Step 3: Overwrite default values as desired

mp.estimator = 'iefc';
%--Whether to perform a model-based (instead of empirical) grid search for the controller
mp.ctrl.flagUseModel = false; 

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;


%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;
mp.runLabel = sprintf('Series%04d_Trial%04d', mp.SeriesNum, mp.TrialNum);

% %--Use just 1 wavelength for initial debugging/testing of code
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
% mp.Nitr = 3; %--Number of wavefront control iterations

% mp.dm1.fourier_spacing = 1.0; %0.5; % Center-to-center spacing between Fourier modes in the focal plane. [lambda/D]
% mp.dm1.fourier_gridType = 'hex';  % Options: 'hex' or 'square'. 'hex' has a denser packing
% xiMin = mp.Fend.corr.Rin-1;
% clocking = mp.Fend.clockAngDeg - 90; % -90 for 'bottom' dark hole.
% [mp.dm1.fourier_basis_xis , mp.dm1.fourier_basis_etas] = falco_choose_fourier_locations_polar(...
%     mp.dm1.Nact/2, mp.dm1.fourier_spacing, mp.dm1.fourier_gridType, xiMin, mp.Fend.corr.Rout+1, mp.Fend.corr.ang, clocking, xiMin);


% % % % PLANNED SEARCH EFC DEFAULTS
% mp.controller = 'plannedEFC';
% mp.ctrl.dmfacVec = 1;
mp.ctrl.log10regVec = -8:1:-4; %--log10 of the regularization exponents (often called Beta values)
% mp.ctrl.log10regVec = -10:1:-3; 

% %--CONTROL SCHEDULE. Columns of mp.ctrl.sched_mat are: 
%     % Column 1: # of iterations, 
%     % Column 2: log10(regularization), 
%     % Column 3: which DMs to use (12, 128, 129, or 1289) for control
%     % Column 4: flag (0 = false, 1 = true), whether to re-linearize
%     %   at that iteration.
%     % Column 5: flag (0 = false, 1 = true), whether to perform an
%     %   EFC parameter grid search to find the set giving the best
%     %   contrast .
%     % The imaginary part of the log10(regularization) in column 2 is
%     %  replaced for that iteration with the optimal log10(regularization)
%     % A row starting with [0, 0, 0, 1...] is for relinearizing only at that time
% mp.ctrl.sched_mat = [...
%     repmat([1, -3, 1, 0, 0], [4, 1]);...
%     ];
% [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);


%--Set path and filename for saved Jacobians
% mp.path.jac = % Define path to saved out Jacobians. Default if left empty is 'falco-matlab/data/jac/'
mp.jac.fn = 'jac_iefc_2.mat'; %'jac_iefc_test.mat'; % Name of the Jacobian file to save or that is already saved. The path to this file is set by mp.path.jac.
mp.relinItrVec = [1, 15, 25, 35]; %[];%1; %[];  %--Correction iterations at which to re-compute the Jacobian. Make an empty vector to load mp.jac.fn
% mp.relinItrVec = []; %[];%1; %[];  %--Correction iterations at which to re-compute the Jacobian. Make an empty vector to load mp.jac.fn

% Use the regular Lyot stop without pinholes for SCC
% mp.P4.compact.mask = mp.P4.compact.maskWithoutPinhole;
% mp.P4.full.mask = mp.P4.full.maskWithoutPinhole;

%% Step 4: Flesh out the rest of the variables

[mp, out] = falco_flesh_out_workspace(mp);

%% Define the probes to use 

dmVsin = sum(mp.dm1.basisCube(:, :, 1:mp.dm1.NbasisModes/2), 3);
figure(31); imagesc(dmVsin); axis xy equal tight; colorbar; colormap parula;
dmVcos = sum(mp.dm1.basisCube(:, :, mp.dm1.NbasisModes/2+1:end), 3);
figure(32); imagesc(dmVcos); axis xy equal tight; colorbar; colormap parula;

mp.iefc.probeCube = zeros(mp.dm1.Nact, mp.dm1.Nact, 2);
mp.iefc.probeCube(:, :, 1) = dmVsin;
mp.iefc.probeCube(:, :, 2) = dmVcos;

mp.iefc.modeCoef = 1e-3 /2; %--Gain coefficient to apply to the normalized DM basis sets for the empirical SCC calibration.
mp.iefc.probeCoef = 2e-4 /2; %--Gain coefficient to apply to the stored probe commands used for IEFC state estimation.
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



%% Step 6: Perform the Wavefront Sensing and Control in FALCO

[mp, out] = falco_wfsc_loop(mp, out);


return
%% Code for doing IEFC WFSC standalone, outside of FALCO



probeCoef0 = 0.01;
dVprobe1 = probeCoef0 * dmVsin;
dVprobe2 = probeCoef0 * dmVcos;

%% Compute response matrix

whichDM = 1;
Nbasis = mp.dm1.NbasisModes;
% jac = zeros(mp.Fend.corr.Npix, Nbasis, mp.jac.Nmode);
jac = zeros(2*mp.Fend.corr.Npix, Nbasis, mp.jac.Nmode);
mp.dm1.V = zeros(mp.dm1.Nact, mp.dm1.Nact);
V0 = mp.dm1.V;
dVcoef = 1e-2; %0.1;
% fn_jac = [mp.path.jac filesep 'iefc_jac_test_fourier_halfLamD.mat'];
fn_jac = [mp.path.jac filesep 'iefc_jac_test_fourier_1em2.mat'];
% fn_jac = [mp.path.jac filesep 'iefc_jac_test_fourier.mat'];
% fn_jac = [mp.path.jac filesep 'iefc_jac_test_fourier_sin_only.mat'];


for iBasis = 1:Nbasis

    fprintf('Empirically obtaining DM response matrix for basis mode: %d/%d\n', iBasis, Nbasis);

    dVmode = dVcoef * mp.dm1.basisCube(:, :, iBasis);
    
    mp.dm1.V = V0 + dVmode;
    evPlus1 = falco_est_delta_intensity(mp, whichDM, dVprobe1,iBasis);
    evPlus2 = falco_est_delta_intensity(mp, whichDM, dVprobe2,iBasis);

    mp.dm1.V = V0 - dVmode;
    evMinus1 = falco_est_delta_intensity(mp, whichDM, dVprobe1,iBasis);
    evMinus2 = falco_est_delta_intensity(mp, whichDM, dVprobe2,iBasis);
    
    mp.dm1.V = V0; % reset
    
    for iMode = 1:mp.jac.Nmode
        jac(1:mp.Fend.corr.Npix, iBasis, iMode) = (evPlus1.Eest(:, iMode) - evMinus1.Eest(:, iMode))/(2*dVcoef);
        jac(mp.Fend.corr.Npix+1:end, iBasis, iMode) = (evPlus2.Eest(:, iMode) - evMinus2.Eest(:, iMode))/(2*dVcoef);
    end

%     % Single Probe
%     mp.dm1.V = V0 + dVmode;
%     evPlus1 = falco_est_delta_intensity(mp, whichDM, dVprobe1);
%     mp.dm1.V = V0 - dVmode;
%     evMinus1 = falco_est_delta_intensity(mp, whichDM, dVprobe1);
%     mp.dm1.V = V0; % reset
%     for iMode = 1:mp.jac.Nmode
%         jac(1:mp.Fend.corr.Npix, iBasis, iMode) = (evPlus1.Eest(:, iMode) - evMinus1.Eest(:, iMode))/(2*dVcoef);
%     end


    
%     E2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
%     E2D(mp.Fend.corr.maskBool) = Ediff;
    figure(30); imagesc(dVmode); axis xy equal tight; colorbar; 
    title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
    drawnow;
    
    E2D = zeros(mp.Fend.Nxi, mp.Fend.Neta);
    E2D(mp.Fend.corr.maskBool) = jac(1:mp.Fend.corr.Npix, iBasis, iMode);
    figure(31); imagesc(log10(abs(E2D))); axis xy equal tight; colorbar; 
    title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
    drawnow;
%     figure(32); imagesc(angle(E2D)); axis xy equal tight; colorbar; colormap hsv;
%     title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
%     drawnow;
end
% 
save(fn_jac, 'jac');

whos jac

%%

% return

%% Check the impact of the probes in the focal plane

p1pad = pad_crop(dVprobe1, 200);
p2pad = pad_crop(dVprobe2, 200);

e1 = fftshift(fft2(fftshift(p1pad)));
figure(41); imagesc(abs(e1)); axis xy equal tight; colorbar; 
figure(42); imagesc(angle(e1)); axis xy equal tight; colorbar; 

e2 = fftshift(fft2(fftshift(p2pad)));
figure(43); imagesc(abs(e2)); axis xy equal tight; colorbar; 
figure(44); imagesc(angle(e2)); axis xy equal tight; colorbar; 


%% Perform Control

load(fn_jac, 'jac')

mp.dm1.V = zeros(mp.dm1.Nact, mp.dm1.Nact);

GG = real(jac' * jac);
normFac = max(diag(GG));

log10regVec = -6:1;%-6:-1;
niVec = zeros(size(log10regVec));


Im = falco_get_summed_image(mp);
ni0 = mean(Im(mp.Fend.corr.maskBool));
fprintf('NI at iteration %d is %.3e\n', 0, ni0);

gainFac = 1;

for Itr = 1:10
    
    probeCoef = 1e-3; %1e-3;%3e-3; %1e-2;%1e-3; %1e-5; %probeCoef0 * 0.9^Itr; %0.75;
    dVprobe1 = probeCoef*dmVsin;
    dVprobe2 = probeCoef*dmVcos;

    whichDM = 1;
    V0 = mp.dm1.V;
    duCube = zeros(mp.dm1.Nact, mp.dm1.Nact, length(log10regVec));

%     % Single probe
%     evProbe1 = falco_est_delta_intensity(mp, whichDM, dVprobe1);
%     DeltaI = evProbe1.Eest;

    evProbe1 = falco_est_delta_intensity(mp, whichDM, dVprobe1);
    evProbe2 = falco_est_delta_intensity(mp, whichDM, dVprobe2);
    DeltaI = ones(2*mp.Fend.corr.Npix, 1);
    DeltaI(1:mp.Fend.corr.Npix) = evProbe1.Eest;
    DeltaI(mp.Fend.corr.Npix+1:end) = evProbe2.Eest;
    
    DeltaI = (probeCoef0/probeCoef) * DeltaI;
    

    for ii = 1:length(log10regVec)
        log10reg = log10regVec(ii);

        duVec = -(GG + normFac*10^(log10reg)*eye(size(GG)))\real(jac.'*DeltaI);
        % https://www.mathworks.com/matlabcentral/answers/400137-computing-a-weighted-sum-of-matrices
        du2D = sum(mp.dm1.basisCube .* reshape(duVec, 1, 1, []), 3);
        duCube(:, :, ii) = du2D;
        figure(41); imagesc(du2D); axis xy equal tight; colorbar; colormap parula; drawnow;

        mp.dm1.V = V0 + du2D;
        Im = falco_get_summed_image(mp);
        figure(42); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, log10(Im), [-9, -3]); axis xy equal tight; colorbar; colormap parula;
        set(gca, 'Fontsize', 20);
        set(gcf, 'Color', 'w');
        drawnow; 
%         pause(0.5);

        niVec(ii) = mean(Im(mp.Fend.corr.maskBool));

    end
    
    figure(43); semilogy(log10regVec, niVec, 'Linewidth', 2);
    set(gca, 'Fontsize', 20);
    set(gcf, 'Color', 'w');
    drawnow;

    
    % Choose best
    [minVal, indexMin] = min(niVec);
    fprintf('NI at iteration %d is %.3e\n', Itr, minVal);
    mp.dm1.V = V0 + duCube(:, :, indexMin);
    
    % Updated image
    Im = falco_get_summed_image(mp);
    figure(44); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, log10(Im), [-9, -3]); axis xy equal tight; colorbar; colormap parula;
    set(gca, 'Fontsize', 20);
    set(gcf, 'Color', 'w');
    drawnow; 

end
