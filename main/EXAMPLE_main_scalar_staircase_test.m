% Copyright 2018-2021 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

clear all
tic;

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command:  addpath(full_path_to_proper); savepath;

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '.'; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
mp.path.ws = '.'; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false

full_path_to_falco = '/Users/niyatid/falco-matlab/';
addpath(genpath(full_path_to_falco));

%% Step 2: Load default model parameters

bws = [0.01,0.01,0.05,0.1,0.15,0.2];
nsbps = [1,3,3,5,7,9];
zernords = [2,3,4,5,6,7,8];
% RMSs = [0.01,0.1,0.5,1,5,10,15,20]; %in nm
RMSs = [0.01,0.1,1,10,50,100,200,300,400,500,550]; %in nm
vals = [];

for index = 1:length(RMSs) %length(bws)
    clearvars -except vals bws index nsbps zernords RMSs
    
    mp.fracBW = bws(6);%index);       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
    mp.Nsbp = nsbps(6);%index);            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
    mp.P1.full.Nbeam = 300; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.P1.compact.Nbeam = 300; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic

    EXAMPLE_defaults_VC_simple
%     EXAMPLE_defaults_SVC_chromatic
    
    mp.F3.full.res = 8; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m
    mp.F3.inVal = 10; % radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m
    mp.F3.outVal = 17;% radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m
    
    mp.F3.VortexCharge = 8;
    mp.F3.NstepStaircase = 8;

    mp.F3.phaseMaskType = 'frenchwrapped';
    mp.F3.VortexSpotDiam = 0.1;

    [mp, out] = falco_flesh_out_workspace(mp);

    %% Calculate image 


%     tic;
%     im = falco_get_summed_image(mp);
%     toc; 
    
%     rawcontrast = mean(im(mp.Fend.score.mask))
%     val = rawcontrast;
%     vals = [vals val];
    
    
    %For Zernike Analysis
    mp.eval.indsZnoll = zernords; %which Zernikes (tip tilt, etc)
    mp.eval.Rsens = [2,4]; %radii to evaluate over
    mp.full.ZrmsVal = RMSs(index)*1E-9;
    sensout = falco_get_Zernike_sensitivities(mp);
    val = sensout;
    
    
    
    vals = [vals val];
end

% %% Plots
% 
% %%-- plot FPM 
% phaseScaleFac = 1;
% pixPerLamD = mp.F3.full.res;
% inputs.type = mp.F3.phaseMaskType;
% inputs.N = ceil_even(pixPerLamD*mp.P1.full.Nbeam);
% inputs.charge = mp.F3.VortexCharge;
% inputs.phaseScaleFac = phaseScaleFac;
% inputs.clocking = mp.F3.clocking;
% inputs.Nsteps = mp.F3.NstepStaircase;
% fpm = falco_gen_azimuthal_phase_mask(inputs); clear inputs;
% 
% figure(1);
% imagesc(abs(fpm));
% colorbar; 
% colormap(gray); 
% 
% figure(2);
% imagesc(angle(fpm));
% colorbar; 
% colormap(hsv);
% caxis([-pi pi]) 
% 
% %%-- plot image 
% 
% figure(3);
% imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(im));
% colorbar; 
% caxis([-12 -5])


% rawcontrast = mean(im(mp.Fend.score.mask))

%% Save data
toc;

%test plot the data being saved
% figure(3)
% bws(1) = 0;
% xaxis = bws;
% plot(xaxis,vals,'Color',[0 0.5 0.8],'LineWidth',2)
% xlabel('Bandwidth');
% ylabel('Raw Contrast');
% title('Raw Contrasts Bandwidth Dependence for SVC')
% set(gca, 'YScale', 'log')

% save classical8contrasts.mat vals bws
save bbfrench8zernikes.mat vals bws zernords RMSs