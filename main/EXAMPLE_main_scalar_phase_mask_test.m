% Copyright 2018-2021 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

clear
tic;

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command:  addpath(full_path_to_proper); savepath;

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = '.'; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = '.'; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false

%% Step 2: Load default model parameters

    
    bws = [0.01,0.01,0.05,0.1,0.15,0.2];
    nsbps = [1,3,3,5,7,9];
    vals = [];
    nulldepths = [];
    peaks =[];
    radii = [0.5 0.65 0.8 0.95 1.1 1.25 1.4 1.55];
    imcube = [];

for index = 1:1  %length(radii) %6 %
    clearvars -except vals bws index nsbps nulldepths peaks radii imcube

    mp.fracBW = bws(1);       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
    mp.Nsbp = nsbps(1);            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
    mp.P1.full.Nbeam = 600; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.P1.compact.Nbeam = 600; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic

    %     EXAMPLE_defaults_VC_simple
    EXAMPLE_defaults_SVC_chromatic
    mp.flagPlot = true;

    mp.Fend.res = 10;
    mp.F3.compact.res = 4;
    mp.F3.full.res = 8; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m
    mp.F3.inVal = 10; % radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m
    mp.F3.outVal = 17;% radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m

    mp.F3.VortexCharge = 6;
%     mp.F3.Nsteps = 6;

    mp.F3.phaseMaskType = 'sawtooth';%plaintest';
    mp.F3.dimpleFlag = false;
    mp.F3.VortexSpotDiamVec = 1.06;%radii(index);%[1.41 1.03];%[1.03 1.41];%1.06;%
    mp.F3.VortexSpotAmpVec = 1;%[1 1];%0;%
    mp.F3.VortexSpotPhaseVec = -0.5;%[0.92 -0.45];%0.5;%

    [mp, out] = falco_flesh_out_workspace(mp);

        %% Calculate image 
% 
% 
%     %     tic;
% %         mp.flagLambdaOffset = true;
% 
%         %For null depth plot
%         %run 20% BW sim once
%         nulldepths = [];
%         phasescalefacs = mp.F3.phaseScaleFac;
%         for iSubband = 1:mp.Nsbp 
%             im3 = falco_get_sbp_image(mp, iSubband);
%             peakheight = max(im3,[],'all');
%             peaks = [peaks peakheight];
%             nulldepths = [nulldepths mean(im3(mp.Fend.score.mask))];
%         end
%     %     toc; 
% % 
        im = falco_get_summed_image(mp);
        imcube(:,:,index) = im;
        %%
%         figure(); imagesc(im/2^16);
%         axis tight; axis equal;
%         set(gca,'ColorScale','log','FontSize',12);
%         h = colorbar; caxis([-7 0]);
        
        xisDL = mp.Fend.xisDL;
        etasDL = mp.Fend.etasDL;

        figure();
%         imagesc(xisDL,etasDL,im);
        imagesc(im/2^16);
        axis image; set(gca,'ydir','normal')
        % caxis([1e-8,5e-5]);
        title("Final Focal Plane Image for Radii = " +radii(index));
        set(gca,'tickdir','out')
        set(gcf,'Color','w');
        colorbar;
        set(gca,'ColorScale','log')
        
        
%         scoremask = im(mp.Fend.score.mask);
%         rawcontrast = mean(im(mp.Fend.score.mask))
%         val = rawcontrast;
%         vals = [vals val];
end

% %% Plots
% % 
% % %%-- plot FPM 
% phaseScaleFac = 1;
% pixPerLamD = mp.F3.full.res;
% inputs.type = mp.F3.phaseMaskType;
% inputs.N = ceil_even(pixPerLamD*mp.P1.full.Nbeam);
% inputs.charge = mp.F3.VortexCharge;
% inputs.phaseScaleFac = phaseScaleFac;
% % inputs.clocking = mp.F3.clocking;
% inputs.Nsteps = mp.F3.Nsteps;
% fpm = falco_gen_azimuthal_phase_mask(inputs); clear inputs;
% % 
% % figure(1);
% % imagesc(abs(fpm));
% % colorbar; 
% % colormap(gray); 
% % 
% figure(3);
% imagesc(angle(fpm));
% colorbar; 
% colormap(hsv);
% caxis([-pi pi]) 
% axis tight; axis equal;
% 
% %%-- plot image 
% 
% figure(3);
% imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(im));
% colorbar; 
% caxis([-12 -5])


% rawcontrast = mean(im(mp.Fend.score.mask))

%% Debugging

modvar = ModelVariables;
mp.debug= true;
[Eout, sDebug] = model_compact(mp,modvar);
figure(6);
imagesc(abs(sDebug.EP4_after_mask));
axis equal; axis tight;


%% Save data
% % toc;

% figure(4)
% bws(1) = 0;
% xaxis = bws;
% xaxis = phasescalefacs;
% plot(xaxis,peaks,'Color',[0 0.5 0.8],'LineWidth',2)
% xlabel('Bandwidth');
% ylabel('Raw Contrast');
% title('Vortex Phase Mask Alone')
% set(gca, 'YScale', 'log')

% save monochromaticroddier.mat vals radii im mp %phasescalefacs mp %imcube mp %nulldepths phasescalefacs peaks im3 mp %