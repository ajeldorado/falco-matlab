% Copyright 2018-2021 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform simple test of SVCs across bandiwdth
% NO DARK HOLE DIGGING- NO WFSC LOOP

% Can test designs with varying roddier radii/phase dimples.

% REVISION HISTORY:
% --------------
% Created on 2022-06-21 by Niyati Desai.


% clear
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
    nsbps = [1,3,3,5,7,9,11];
    vals = [];
    nulldepths = [];
    peaks =[];
    % radii = [0.5 0.65 0.8 0.95 1.1 1.25 1.4 1.55];
    % radiilong = linspace(0,0.7,7); 
    radii = [0.3]; %linspace(0.01,0.7,7); %linspace(0,0.2,7);
    phaselist = [0 0.2 0.4 0.6];
    imcube = [];
    
    
    seriesnum = 2;
    trialnum = 5;
    flagSaveOutput = false;


for index = 1:length(radii) %6 %
    clearvars -except flagSaveOutput spotvals holevals seriesnum trialnum vals radii bws index nsbps nulldepths peaks phaselist imcube

    %--Record Keeping
    mp.SeriesNum = seriesnum;
    mp.TrialNum = trialnum;

    mp.fracBW = bws(1);       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
    mp.Nsbp = nsbps(1);            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
    mp.P1.full.Nbeam = 300; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.P1.compact.Nbeam = 300; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic

%         EXAMPLE_defaults_VC_simple
    EXAMPLE_defaults_SVC_chromatic
    mp.flagPlot = false;

    mp.Fend.res = 35.55;
    mp.F3.compact.res = 16;
    mp.F3.full.res = 16; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m
    mp.F3.inVal = 1; % radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m
    mp.F3.outVal = 2;% radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m

    mp.F3.phaseMaskType = 'roddier'; %'sawtooth'; %'roddier';
    mp.F3.VortexCharge = 6;
    mp.F3.roddierradius = 0.621; %[lambda/D]
    mp.F3.roddierphase = 0.5; %phaselist(index);
    
    if strcmp(mp.F3.phaseMaskType,'roddier')
        mp.F3.flagDimple = true; %true;
    else
        mp.F3.flagDimple = false;
    end
    lams = 1e6.*mp.lambda0*linspace(1-mp.fracBW/2, 1+mp.fracBW/2, mp.Nsbp);
    
    mp.flagAchromat = false;
    %Uncomment this line if testing dual-layer achromat
    if mp.flagAchromat
        addpath(genpath('/Users/ndesai/Documents/falco-matlab/lib-external/materials'));
        mp.F3.material1 = "LiNbO3";
        mp.F3.material2 = "TiO2";
        [mp.F3.thickness,mp.F3.phaseScaleFac] = optimizeachromat(mp.F3.material1,mp.F3.material2,lams);
    end
    mp.flagVVC = false; %false;

    if mp.flagVVC
        fprintf('achromatic VVC')
        % mp.sbp_weights = ones(mp.Nsbp,1);
        % mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBW/2, 1+mp.fracBW/2, mp.Nsbp);
        % mp.sbp_weights(1) = 1/2; %--Give end sub-bands half weighting
        % mp.sbp_weights(end) = 1/2; %--Give end sub-bands half weighting
        % mp.F3.phaseScaleFacLambdas = ones(1, mp.Nsbp) * mp.lambda0;
        % mp.F3.phaseScaleFac = ones(1, mp.Nsbp);
        mp.F3.phaseScaleFac = 1;
    end

    %Mask Errors
    
    mp.F3.full.flagErrors = true;
    mp.F3.holeradius = radii(index);
    
    % r = mp.F3.holeradius*resErrorMap; %20; 
    
    % Create array with central circular region
    Error_res_full = floor(mp.F3.full.res*mp.P1.full.Nbeam/(2*mp.F3.outVal));
    inputs.pixresFPM = Error_res_full; %--pixels per lambda/D
    beamDiamPix = ceil_even(mp.F3.full.res*mp.P1.full.Nbeam); %mp.P1.full.Nbeam;
    inputs.rhoInner = mp.F3.holeradius; % radius of inner FPM amplitude spot (in lambda_c/D)
    inputs.rhoOuter = inf; % radius of outer opaque FPM ring (in lambda_c/D)
    inputs.FPMampFac = 0; % amplitude transmission of inner FPM spot
    inputs.centering = 'pixel';
    EerrorMap = falco_gen_annular_FPM(inputs);
    EerrorMap = pad_crop(EerrorMap, beamDiamPix, 'extrapval', 1);

    mp.F3.full.EabArray = EerrorMap .* reshape(mp.F3.phaseScaleFac, 1, 1, []);
    mp.F3.full.EabRes = Error_res_full;

    % Error_res_compact = floor(mp.F3.compact.res*mp.P1.compact.Nbeam/(2*mp.F3.outVal));
    % inputs.pixresFPM = Error_res_compact;
    % beamDiamPix_compact = ceil_even(mp.F3.compact.res*mp.P1.compact.Nbeam);
    % EerrorMap_compact = falco_gen_annular_FPM(inputs);
    % EerrorMap_compact = pad_crop(EerrorMap_compact, beamDiamPix_compact, 'extrapval', 1);
    % 
    % mp.F3.compact.EabArray = EerrorMap_compact .* reshape(mp.F3.phaseScaleFac, 1, 1, []);
    % mp.F3.compact.EabRes = Error_res_compact;

    figure(12);
    imagesc(abs(EerrorMap));
    axis xy equal tight; colorbar;
    colormap gray;
    
    figure(13);
    imagesc(angle(EerrorMap));
    axis xy equal tight; colorbar;
    colormap hsv;
    
    [mp, out] = falco_flesh_out_workspace(mp);

        %% Calculate image 


    %     tic;
%         mp.flagLambdaOffset = true;

        %For null depth plot
        %run 20% BW sim once
        nulldepths = [];
        phasescalefacs = mp.F3.phaseScaleFac;
        for iSubband = 1:mp.Nsbp 
            im3 = falco_get_sbp_image(mp, iSubband);
%             peakheight = max(im3,[],'all');
%             peaks = [peaks peakheight];
            nulldepths = [nulldepths mean(im3(mp.Fend.score.mask))];

        end
%     %     toc; 
% % 
        im = falco_get_summed_image(mp);
%         imcube(:,:,index) = im;
 
%         figure();
% %         imagesc(xisDL,etasDL,im);
%         imagesc(im/2^16);
%         axis image; set(gca,'ydir','normal')
%         % caxis([1e-8,5e-5]);
%         title("Final Focal Plane Image");
%         set(gca,'tickdir','out')
%         set(gcf,'Color','w');
%         colorbar;
%         set(gca,'ColorScale','log')
%         
% IF ITERATING OVER SOME VARIABLE (other than broadband)         
        scoremask = im(mp.Fend.score.mask);
        rawcontrast = mean(im(mp.Fend.score.mask))
        val = rawcontrast;
        vals = [vals val];

end

%% Plots

% 
%%-- plot image 
% 
% figure(5);
% imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(im));
% colorbar; 
% caxis([-15 -5]);
% axis equal; axis tight;
% title(['Final Focal Plane Image'])
% set(gca,'tickdir','out')
% set(gcf,'Color','w');
% colorbar;
% set(gca,'ColorScale','log')
%%

% rawcontrast = mean(im(mp.Fend.score.mask))
% thicknesses = mp.F3.thickness

%% Debugging

% modvar = ModelVariables;
% mp.debug= true;
% [Eout, sDebug] = model_compact(mp,modvar);
% figure(6);
% imagesc(abs(sDebug.EP4_after_mask));
% axis equal; axis tight;


%% Save data
spotvals = vals;
% holevals = vals;
% spotvalslong = vals;
fnumber = 83; %32;
% bws(1) = 0;
% xaxis = radiilong*83*760e-3;
% 
figure(49)
% xaxis = mp.sbp_centers; %phasescalefacs;
plot(radii*fnumber*mp.lambda0*1e6,spotvals,'Color',[0 0.5 0.8],'LineWidth',2)
% hold on;
% plot(radii*fnumber*mp.lambda0*1e6,holevals,'Color',[0.2 0.1 0.8],'LineWidth',2)
% hold on;
% plot(xaxis,spotvalslong,'Color',[0 0.5 0.8],'LineWidth',2)

xlabel('Spot size [radius in microns]'); %'\lambda');%/\lambda_0');
ylabel('Raw Contrast from 3-10 $\lambda$/D');
title('1% BW Scalar Vortex WITHOUT BABINET'); %'Diamond+TiO2 sawtooth')
set(gca, 'YScale', 'log')
legend('SVC with Opaque Spot')
% legend('Spot','Hole')
grid on;

% out.fnumber = fnumber;
% out.wavs = mp.sbp_centers;
% out.radii = radii;
% % out.nulldepths = nulldepths;
% out.mp = mp;
% % out.spotvals = spotvals;
% out.holevals = holevals;
% % out.spotvalslong = spotvalslong;
% out.finalim = im;

dataDir = '/Users/ndesai/Documents/Sim_Data/';
seriesX = mp.SeriesNum;
trialY = mp.TrialNum;
%%
if flagSaveOutput
    addpath '/Users/ndesai/Documents/Sim_Data';
    saveTrialData(dataDir, seriesX, trialY, out)
end
% save filenamehere.mat vals radii im mp %phasescalefacs mp %imcube mp %nulldepths phasescalefacs peaks im3 mp %