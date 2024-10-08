% Copyright 2022, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a loop of SVC simulation runs WITH WFSC!!
%
% Can be used for Contrast measurements
% for varying bandwidths, resolutions, and radii params
%
% REVISION HISTORY:
% --------------
% Created on 2022-02-23 by Niyati Desai.
% Added roddier capabilities 2022-011-23
% Added metasurface capabilities 2023-05-24
 
clear
close all
clear all
 

bws = [0.01,0.01,0.05,0.1,0.15,0.8/3.8];
nsbps = [1,3,3,5,7,9];
vals = [];
res = [100,200,300,400,500,600,700]; %found 600 is sufficient

radiiList = linspace(0.2,1.8,17);
phaseList = linspace(0,1,11);
tic; 
 
for index = 1:1 %length(RMSs) %length(res)
    toc;
    clearvars -except vals bws index nsbps res 
    mp.use_lastJacStruc = false;
    
    %% Step 1: Define Necessary Paths on Your Computer System
    slowpoke = false;
    
    if ~slowpoke
        % %--Library locations. FALCO and PROPER are required. CVX is optional.
        % mp.path.falco = '/Users/ndesai/Documents/falco-matlab/';  %--Location of FALCO
        % 
        % %%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
        % mp.path.config = '/Users/ndesai/Documents/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
        % mp.path.ws = '/Users/ndesai/Documents/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
        % mp.path.mask = '/Users/ndesai/Documents/falco-matlab/lib/masks/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
        % mp.path.ws_inprogress = mp.path.ws;
    else
    % for slowpoke machine
 
        % % %--Library locations. FALCO and PROPER are required. CVX is optional.
        mp.path.falco = 'C:\Users\jdllop\Documents\GitHub\falco-matlab';%'~/Repos/falco-matlab/';  %--Location of FALCO
 
        % %%--Output Data Directories ( Comment these lines out to use defaults within falco-matlab/data/ directory.)
        mp.path.config = 'C:\Users\jdllop\Documents\GitHub\falco-matlab\data\brief';%'~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
        mp.path.ws = 'C:\Users\jdllop\Documents\GitHub\falco-matlab\data\ws';%'~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
        
        %%--Add to the MATLAB Path
        addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
    
    end
 
    %% Step 2: Load default model parameters
 
    disp(index);
    mp.fracBW = bws(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.Nsbp = nsbps(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.P1.full.Nbeam = 300; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.P1.compact.Nbeam = 300; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    
    EXAMPLE_defaults_SVC_chromatic
%     EXAMPLE_defaults_VC_simple
    
    % Specify SVC type: Valid options are 'vortex', 'cos', 'sectors', 'staircase', 'sawtooth', 'wrapped8', 'wrapped6', 'dzpm', 'roddier', and 'just_dimple' 
    mp.F3.phaseMaskType = 'roddier';
    mp.F3.VortexCharge = 6;
    mp.F3.NstepStaircase = 6;
    mp.F3.roddierradius = 0.53; %[lambda/D]
    mp.F3.roddierphase = 0.5;
    mp.F3.flagDimple = true;
    
 
 
    %% Step 3: Overwrite default values as desired
    mp.Fend.res = 35.5;
    mp.F3.compact.res = 6; %4;
    mp.F3.full.res = 16; %8; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m
    mp.F3.inVal = 0; %10; % radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m
    mp.F3.outVal = 5; % 17;% radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m
 
    %%--Special Computational Settings
    mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
    mp.flagPlot = true;
 
    %--Record Keeping
    mp.SeriesNum = 1;
    mp.TrialNum = 1;
 
    mp.Nwpsbp = 1;          %--Number of wavelengths to be used to approximate an image in each sub-bandpass
    mp.Nitr = 2; %--Number of wavefront control iterations
 
    mp.flagVVC = false;

    % progress plothe 
    mp.plot.contrast_clim_min = -15;
    mp.plot.contrast_clim_max = -3; 

    if mp.flagVVC
        fprintf('achromatic VVC')
        mp.sbp_weights = ones(mp.Nsbp,1);
        mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBW/2, 1+mp.fracBW/2, mp.Nsbp);
        mp.sbp_weights(1) = 1/2; %--Give end sub-bands half weighting
        mp.sbp_weights(end) = 1/2; %--Give end sub-bands half weighting
        mp.F3.phaseScaleFacLambdas = ones(1, mp.Nsbp) * mp.lambda0;
        mp.F3.phaseScaleFac = ones(1, mp.Nsbp);
    end
    
    %% Step 4: Generate the label associated with this trial
 
    mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
        mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
        '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
        '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
        '_',mp.controller];
 
 
    %% Step 5: Perform the Wavefront Sensing and Control
 
    [mp, out] = falco_flesh_out_workspace(mp);
    
    
    %For Contrasts with WFSC
    [mp, out] = falco_wfsc_loop(mp, out);
    val = out.InormHist(end);
    
    
    %For Contrast w/o WFSC (0 wavefront control iterations)
%     outSingle = falco_eval_without_control(mp);
%     val = outSingle.InormHist(1);
      
      vals = [vals val];
      toc;
 
 
end

 
figure()
ax = gca;
 
plot(0:length(out.InormHist)-1,out.InormHist,'-^','color',[0.4660 0.6740 0.1880],'LineWidth',2);
grid on;
title('20% Broadband Convergence Rate','fontsize',19,'Color','w')
% xlabel('Iteration','fontsize',19,'Color','w')
% ylabel('Mean Normalized Intensity','fontsize',19,'Color','w')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
xlabel('Iteration','fontsize',19,'Color','k')
ylabel('Mean Normalized Intensity','fontsize',19,'Color','k')
set(gca, 'LineWidth',1);
set(gca, 'YScale', 'log','fontsize',19)
% set(gca,'Color','k')%color for the plot area
% set(gca,'XColor',[1 1 1]); % Set RGB value to what you want
% set(gca,'YColor',[1 1 1]); % Set RGB value to what you want
% set(gcf,'color','k');
set(gcf,'color','w');
grid on;
box off;
 
%%
 
nulldepths = [];
phasescalefacs = mp.F3.phaseScaleFac;
for iSubband = 1:mp.Nsbp 
    im3 = falco_get_sbp_image(mp, iSubband);
%             peakheight = max(im3,[],'all');
%             peaks = [peaks peakheight];
    nulldepths = [nulldepths mean(im3(mp.Fend.score.mask))]
end
 
% 
im = falco_get_summed_image(mp);
 
%         
xisDL = mp.Fend.xisDL;
etasDL = mp.Fend.etasDL;
%% 
figure();
    imagesc(xisDL,etasDL,im);
    axis xy equal tight;
    caxis([1e-11,5e-6]);
    title("Final Focal Plane Image for Phase");
    set(gca,'tickdir','out')
    set(gcf,'Color','w');
    colorbar;
%     set(gca,'ColorScale','log')
        
figure(44)
    % bws(1) = 0;
    % xaxis = bws;
    xaxis = phasescalefacs;
    plot(xaxis,nulldepths,'Color',[0 0.5 0.8],'LineWidth',2)
    xlabel('Bandwidth');
    ylabel('Raw Contrast');
    title('Chromatic Performance')
    set(gca, 'YScale', 'log')
    grid on;
 
%%
figure()
imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(im),[-8.5 -3]); 
axis xy equal tight; colorbar; colormap(parula);
%(h_psf); colormap(h_psf,parula);
% xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',20,'FontName','Times','FontWeight','Normal')
set(gcf,'color','w')
%ylabel(ch_psf,'$log_{10}$(NI)','Fontsize',24,'Interpreter','LaTex');
title(sprintf('Stellar PSF: NI = %.2e',out.InormHist(end)),'Fontsize',24);%,'Fontweight','Bold');
 

% save frenchwrappedzernikes.mat vals bws zernords RMSs
