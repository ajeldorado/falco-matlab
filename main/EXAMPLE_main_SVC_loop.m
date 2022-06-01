% Copyright 2022, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform a loop of SVC simulation runs.
%
% Can be used for Contrast measurements or Zernike analysis
% for varying bandwidths, RMSs, and resolutions


% REVISION HISTORY:
% --------------
% Created on 2022-02-23 by Niyati Desai.


clear all;
close all;


RMSs = [0.01,0.1,0.5,1,5,10,15,20]; %in nm
bws = [0.01,0.01,0.05,0.1,0.15,0.2];
nsbps = [1,3,3,5,7,9];
vals = [];
res = [100,200,300,400,500,600,700]; %found 600 is sufficient
zernords = [2,3,4,5,6,7,8];


for index = 1:2;%length(RMSs) %length(res)
    clearvars -except vals bws index nsbps res RMSs zernords
    mp.use_lastJacStruc = false;
    
    %% Step 1: Define Necessary Paths on Your Computer System
    slowpoke = false;
    
    if ~slowpoke
        %--Library locations. FALCO and PROPER are required. CVX is optional.
        mp.path.falco = '/Users/niyatid/falco-matlab/';  %--Location of FALCO
        mp.path.proper = '/Users/niyatid/falco-matlab/lib_external/proper/'; %--Location of the MATLAB PROPER library
        % mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

        %%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
        mp.path.config = '/Users/niyatid/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
        mp.path.ws = '/Users/niyatid/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
        mp.path.mask = '/Users/niyatid/falco-matlab/lib/masks/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
        mp.path.ws_inprogress = mp.path.ws;
    else
    %for slowpoke

        % % %--Library locations. FALCO and PROPER are required. CVX is optional.
        mp.path.falco = 'C:\Users\jdllop\Documents\GitHub\falco-matlab';%'~/Repos/falco-matlab/';  %--Location of FALCO
        mp.path.proper = 'C:\Users\jdllop\Documents\GitHub\falco-matlab\proper';%'~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library

        % %%--Output Data Directories ( Comment these lines out to use defaults within falco-matlab/data/ directory.)
        mp.path.config = 'C:\Users\jdllop\Documents\GitHub\falco-matlab\data\brief';%'~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
        mp.path.ws = 'C:\Users\jdllop\Documents\GitHub\falco-matlab\data\ws';%'~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];
    end
 
    %%--Add to the MATLAB Path
    addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
    addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


    %% Step 2: Load default model parameters

    disp(index);
    mp.fracBW = bws(1);%index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.Nsbp = nsbps(1);%index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.P1.full.Nbeam = 600; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.P1.compact.Nbeam = 600; %res(index); %make sure this line is commented out in EXAMPLE_defaults_HCST_SVC_chromatic
    mp.F3.phaseMaskType = 'frenchwrapped';
%     mp.F3.VortexCharge = -8;
%     mp.F3.NstepStaircase = 6;
    EXAMPLE_defaults_HCST_SVC_chromatic


    %% Step 3: Overwrite default values as desired
    
    mp.F3.full.res = 8; % Coarse DFT resolution used in propcustom_mft_PtoFtoP.m
    mp.F3.inVal = 10; % radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m
    mp.F3.outVal = 17;% radius of fine-sampled DFT region in propcustom_mft_PtoFtoP.m
    mp.F3.VortexSpotDiam = 0.1;

    %%--Special Computational Settings
    mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
    mp.flagPlot = true;

    %--Record Keeping
    mp.SeriesNum = 1;
    mp.TrialNum = 1;

    mp.Nwpsbp = 1;          %--Number of wavelengths to be used to approximate an image in each sub-bandpass
    mp.Nitr = 1; %--Number of wavefront control iterations

    %% Step 4: Generate the label associated with this trial

    mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
        mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
        '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
        '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
        '_',mp.controller];


    %% Step 5: Perform the Wavefront Sensing and Control

    [mp, out] = falco_flesh_out_workspace(mp);

    
%     %For Zernike Analysis
%     mp.eval.indsZnoll = zernords; %which Zernikes (tip tilt, etc)
%     mp.eval.Rsens = [2,4]; %radii to evaluate over
%     mp.full.ZrmsVal = RMSs(index)*1E-9;
%     sensout = falco_get_Zernike_sensitivities(mp);
%     val = sensout;
    
    
    %For Contrasts with WFSC
%     [mp, out] = falco_wfsc_loop(mp, out);
%     val = out.InormHist(end);
    
    
    %For Contrast w/o WFSC (0 wavefront control iterations)
%     outSingle = falco_eval_without_control(mp);
%     val = outSingle.InormHist(1);
      
      vals = [vals val];


end
% save frenchwrappedzernikes.mat vals bws zernords RMSs
