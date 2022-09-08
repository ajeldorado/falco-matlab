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
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false


%% Step 2: Load default model parameters

EXAMPLE_defaults_SCC


%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

% %--Use just 1 wavelength for initial debugging/testing of code
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
% mp.Nitr = 3; %--Number of wavefront control iterations


%% Step 4: Generate the label associated with this trial

mp.runLabel = sprintf('Series%04d_Trial%04d', mp.SeriesNum, mp.TrialNum);


%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);


%%
return

%% 
mp.P4.full.mask = mp.P4.full.maskWithPinhole;
[mp, out] = falco_flesh_out_workspace(mp);
ev1 = falco_est_scc(mp);

% mp.P4.full.mask = mp.P4.full.maskWithoutPinhole;
% [mp, out] = falco_flesh_out_workspace(mp);
% ev2 = falco_est_perfect_Efield_with_Zernikes(mp);
% 
% 
% E2D1 = zeros(mp.Fend.Neta, mp.Fend.Nxi);
% E2D1(mp.Fend.corr.maskBool) = ev1.Eest;
% 
% E2D2 = zeros(mp.Fend.Neta, mp.Fend.Nxi);
% E2D2(mp.Fend.corr.maskBool) = ev2.Eest;
% 
% figure(11); imagesc(log10(abs(E2D1))); axis xy equal tight; colorbar; drawnow;
% figure(12); imagesc(angle(E2D1)); axis xy equal tight; colorbar; colormap hsv; drawnow;
% figure(13); imagesc(real(E2D1)); axis xy equal tight; colorbar; drawnow;
% figure(14); imagesc(imag(E2D1)); axis xy equal tight; colorbar; drawnow;
% 
% figure(21); imagesc(log10(abs(E2D2))); axis xy equal tight; colorbar; drawnow;
% figure(22); imagesc(angle(E2D2)); axis xy equal tight; colorbar;  colormap hsv; drawnow;
% figure(23); imagesc(real(E2D2)); axis xy equal tight; colorbar; drawnow;
% figure(24); imagesc(imag(E2D2)); axis xy equal tight; colorbar; drawnow;
% 
% % figure(25); imagesc(abs(E2D2)./abs(E2D1)); axis xy equal tight; colorbar; drawnow;


%% Compute response matrix

Nbasis = mp.dm1.NbasisModes;
jac = zeros(mp.Fend.corr.Npix, Nbasis, mp.jac.Nmode);
mp.dm1.V = zeros(mp.dm1.Nact, mp.dm1.Nact);
V0 = mp.dm1.V;
dVcoef = 1;

fn_jac = [mp.path.jac filesep 'scc_jac_test_fourier.mat'];

% for iBasis = 1:Nbasis
% 
%     fprintf('Empirically obtaining DM response matrix for basis mode: %d/%d\n', iBasis, Nbasis);
%     dV = mp.dm1.basisCube(:, :, iBasis);
%     mp.dm1.V = V0 + dV;
%     evTempPos = falco_est_scc(mp);
%     mp.dm1.V = V0 - dV;
%     evTempNeg = falco_est_scc(mp);
% 
%     for iMode = 1:mp.jac.Nmode
%         Ediff = (evTempPos.Eest(:, iMode) - evTempNeg.Eest(:, iMode))/(2*dVcoef);
%         jac(:, iBasis, iMode) = Ediff;
%     end
% 
% %     E2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
% %     E2D(mp.Fend.corr.maskBool) = Ediff;
% %     figure(30); imagesc(dV); axis xy equal tight; colorbar; 
% %     title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
% %     drawnow;
% %     figure(31); imagesc(log10(abs(E2D))); axis xy equal tight; colorbar; 
% %     title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
% %     drawnow;
% %     figure(32); imagesc(angle(E2D)); axis xy equal tight; colorbar; colormap hsv;
% %     title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
% %     drawnow;
% end
% 
% save(fn_jac, 'jac');


% % Nbasis = mp.dm1.NactTotal;
% % jac = zeros(mp.Fend.corr.Npix, Nbasis, mp.jac.Nmode);
% % mp.dm1.V = zeros(mp.dm1.Nact, mp.dm1.Nact);
% % V0 = mp.dm1.V;
% % dVcoef = 1;
% % fn_jac = [mp.path.jac filesep 'scc_jac_test.mat'];
% %
% % for iBasis = 1:Nbasis
% % 
% %     fprintf('Empirically obtaining DM response matrix for basis mode: %d/%d\n', iBasis, Nbasis);
% % 
% %     dV = zeros(mp.dm1.Nact, mp.dm1.Nact);
% %     dV(iBasis) = dVcoef;
% %     mp.dm1.V = V0 + dV;
% %     evTempPos = falco_est_scc(mp);
% %     mp.dm1.V = V0 - dV;
% %     evTempNeg = falco_est_scc(mp);
% % 
% %     for iMode = 1:mp.jac.Nmode
% %         Ediff = (evTempPos.Eest(:, iMode) - evTempNeg.Eest(:, iMode))/(2*dVcoef);
% %         jac(:, iBasis, iMode) = Ediff;
% %     end
% % 
% % %     E2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
% % %     E2D(mp.Fend.corr.maskBool) = Ediff;
% % %     figure(30); imagesc(dV); axis xy equal tight; colorbar; 
% % %     title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
% % %     drawnow;
% % %     figure(31); imagesc(log10(abs(E2D))); axis xy equal tight; colorbar; 
% % %     title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
% % %     drawnow;
% % %     figure(32); imagesc(angle(E2D)); axis xy equal tight; colorbar; colormap hsv;
% % %     title(sprintf('Basis Mode %d/%d', iBasis, Nbasis));
% % %     drawnow;
% % end
% % 
% % save(fn_jac, 'jac');
% %
% % load(fn_jac)
% % 
% % jacSum = squeeze(sum(abs(jac).^2, 1));
% % jac2D = reshape(jacSum, [mp.dm1.Nact, mp.dm1.Nact]);
% % figure(40); imagesc(log10(jac2D)); axis xy equal tight; colorbar; colormap parula;


%%

load(fn_jac, 'jac')

GG = real(jac' * jac);
normFac = max(diag(GG));

log10regVec = -6:-1;%-6:-1;
niVec = zeros(size(log10regVec));



mp.dm1.V = zeros(mp.dm1.Nact, mp.dm1.Nact);

Im = falco_get_summed_image(mp);
ni0 = mean(Im(mp.Fend.corr.maskBool));
fprintf('NI at iteration %d is %.3e\n', 0, ni0);


for Itr = 1:5

    ev = falco_est_scc(mp);
    V0 = mp.dm1.V;
    duCube = zeros(mp.dm1.Nact, mp.dm1.Nact, length(log10regVec));

    for ii = 1:length(log10regVec)
        log10reg = log10regVec(ii);


        duVec = -(GG + normFac*10^(log10reg)*eye(size(GG)))\real(jac'*ev.Eest);
        % https://www.mathworks.com/matlabcentral/answers/400137-computing-a-weighted-sum-of-matrices
        du2D = sum(mp.dm1.basisCube .* reshape(duVec, 1, 1, []), 3);
        duCube(:, :, ii) = du2D;
        figure(41); imagesc(du2D); axis xy equal tight; colorbar; colormap parula; drawnow;

        mp.dm1.V = V0 + du2D;
        Im = falco_get_summed_image(mp);
        figure(42); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, log10(Im), [-8, -3]); axis xy equal tight; colorbar; colormap parula;
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
    
    

    
end

%% Check NI without Pinhole

mp.P4.full.mask = mp.P4.full.maskWithoutPinhole;
[mp, out] = falco_flesh_out_workspace(mp);
% ev2 = falco_est_perfect_Efield_with_Zernikes(mp);

Im = falco_get_summed_image(mp);
niNoPH = mean(Im(mp.Fend.corr.maskBool));
fprintf('NI without pinhole is %.3e\n', niNoPH);



