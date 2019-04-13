% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to:
%  1) Specify the key parameter values for a hybrid Lyot coronagraph.
%  2) Load the rest of the default settings.
%  3) Save out all the input parameters.
%  4) Run a single trial of WFSC using FALCO.

close all
clear all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

%% Step 2: Load default model parameters

EXAMPLE_defaults_WFIRST_HLC_BMC_WFE

%% Properties for BMC Analysis

mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;           %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control

% mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 6;           %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control

mp.d_dm1_dm2 = 0.6; % distance between DM1 and DM2 (meters)

%%--Pupil Masks
mp.P1.compact.Nbeam = 200;%927;
mp.P2.compact.Nbeam = mp.P1.compact.Nbeam ;
mp.P3.compact.Nbeam = mp.P1.compact.Nbeam ;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam ;

mp.P1.full.Nbeam = 927; %--Number of pixel widths across the actual diameter of the beam/aperture (independent of beam centering)
mp.P2.full.Nbeam = mp.P1.full.Nbeam;
mp.P3.full.Nbeam = mp.P1.full.Nbeam;
mp.P4.full.Nbeam = mp.P1.full.Nbeam;

%--WFE maps and stops on DMs:
mp.dm1.inf_fn = 'influence_BMC_2kDM_400micron_res10.fits';
mp.dm2.inf_fn = 'influence_BMC_2kDM_400micron_res10.fits';

mp.flagDMwfe = true;
if(mp.flagDMwfe)
    mp.dm1.wfe = fitsread(sprintf('~/Data/BMC/wfe/bmc50_dm_wfe_%dpix_pupil.fits',mp.P1.full.Nbeam));
    mp.dm2.wfe = mp.dm1.wfe; %rot90(mp.dm1.wfe,2);
end

% mp.dm1.dm_spacing = 400e-6; %--User defined actuator pitch
% mp.dm2.dm_spacing = 400e-6; %--User defined actuator pitch

mp.flagDM1stop = false; %--logical flag whether to include the stop at DM1 or not
mp.flagDM2stop = true; %--logical flag whether to include the stop at DM2 or not
mp.dm2.Dstop = 49*mp.dm2.dm_spacing;  %--diameter of circular stop at DM2 and centered on the beam

mp.P2.D = 46.2937*mp.dm1.dm_spacing; % beam diameter at pupil closest to the DMs  (meters)
mp.P3.D = mp.P2.D;
mp.P4.D = mp.P2.D;

%mp.dm1.Nact = 50; % number of actuators across DM1
%mp.dm2.Nact = 50; % number of actuators across DM2

%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = false;%true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 36;
mp.TrialNum = 1;

%--Control
mp.Nitr = 50; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

%% [OPTIONAL] Start from a previous FALCO trial's DM settings

% fn_prev = sprintf('Series0033_Trial%04d_HLC_WFIRST180718_3DM50_z%s_IWA2.7_OWA10_6lams575nm_BW10_plannedEFC_snippet.mat',13,num2str(mp.d_dm1_dm2));
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Initialization function (copied from falco_wfsc_loop.m)

%--Sort out file paths and save the config file    

%--Add the slash or backslash to the FALCO path if it isn't there.
if( strcmp(mp.path.falco(end),'/')==false || strcmp(mp.path.falco(end),'\')==false )
    mp.path.falco = [mp.path.falco filesep];
end

mp.path.dummy = 1; %--Initialize the folders structure in case it doesn't already exist

%--Store minimal data to re-construct the data from the run: the config files and "out" structure after a trial go here
if(isfield(mp.path,'config')==false)
    mp.path.config = [mp.path.falco filesep 'data' filesep 'config' filesep];     
end

%--Entire final workspace from FALCO gets saved here.
if(isfield(mp.path,'ws')==false)
    mp.path.ws = [mp.path.falco filesep 'data' filesep 'ws' filesep];      
end

%--Save the config file
fn_config = [mp.path.config mp.runLabel,'_config.mat'];
save(fn_config)
fprintf('Saved the config file: \t%s\n',fn_config)

%--Get configuration data from a function file
[mp,out] = falco_init_ws(fn_config);

%%
%--Update the number of elements used per DM
if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); end
if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); end
if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); end
if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); end

%% Part 5: Compute the Jacobian using model_Jacobian

modvar.wpsbpIndex = mp.wi_ref;
modvar.whichSource = 'star'; 

%--Re-initialize the Jacobian arrays to full size
G1=zeros(1,1,mp.jac.Nmode); G2=zeros(1,1,mp.jac.Nmode); G3=zeros(1,1,mp.jac.Nmode); G4=zeros(1,1,mp.jac.Nmode); G5=zeros(1,1,mp.jac.Nmode); G6=zeros(1,1,mp.jac.Nmode); G7=zeros(1,1,mp.jac.Nmode); G8=zeros(1,1,mp.jac.Nmode); G9=zeros(1,1,mp.jac.Nmode); %--Initialize for bookkeeping in cells later 
if(any(mp.dm_ind==1)); G1 = zeros(mp.Fend.corr.Npix,mp.dm1.NactTotal,mp.jac.Nmode); end % control Jacobian for DM1
if(any(mp.dm_ind==2)); G2 = zeros(mp.Fend.corr.Npix,mp.dm2.NactTotal,mp.jac.Nmode); end % control Jacobian for DM2
if(any(mp.dm_ind==8)); G8 = zeros(mp.Fend.corr.Npix,mp.dm8.NactTotal,mp.jac.Nmode); end % control Jacobian for DM8
if(any(mp.dm_ind==9)); G9 = zeros(mp.Fend.corr.Npix,mp.dm9.NactTotal,mp.jac.Nmode); end % control Jacobian for DM9
    
%--Compute the number of total actuators for all DMs used. 

GallCell1 = {squeeze(G1(:,:,1)),squeeze(G2(:,:,1)),squeeze(G3(:,:,1)),squeeze(G4(:,:,1)),squeeze(G5(:,:,1)),squeeze(G6(:,:,1)),squeeze(G7(:,:,1)),squeeze(G8(:,:,1)),squeeze(G9(:,:,1))}; % Create the cell array. Placeholders for non-existent Jacobians to have consistent numbering
NeleAll = 0;
NeleVec = []; %--Vector of total number of used actuators for each used DM
for ii=1:numel(mp.dm_ind)
    dm_index = mp.dm_ind(ii);
    NeleAll = NeleAll + size(GallCell1{dm_index},2);
    NeleVec = [NeleVec; size(GallCell1{dm_index},2) ];
end
clear GallCell1 %--Save RAM

%--Compute the control Jacobians for each DM
jacStruct =  model_Jacobian(mp);
if(any(mp.dm_ind==1)); G1_jac = jacStruct.G1; end
if(any(mp.dm_ind==2)); G2_jac = jacStruct.G2; end
if(any(mp.dm_ind==8)); G8_jac = jacStruct.G8; end
if(any(mp.dm_ind==9)); G9_jac = jacStruct.G9; end
clear jacStruct  %--Save RAM

%% Part 6: Compute the Jacobian using model_compact and differencing

Vfrac = 1e-6; %--Want delta voltage to be tiny for DM1 and DM2 to stay linear

%--Re-initialize the Jacobian arrays to full size
if(any(mp.dm_ind==1)); G1 = zeros(mp.Fend.corr.Npix,mp.dm1.NactTotal,mp.jac.Nmode); end % control Jacobian for DM1
if(any(mp.dm_ind==2)); G2 = zeros(mp.Fend.corr.Npix,mp.dm2.NactTotal,mp.jac.Nmode); end % control Jacobian for DM2
if(any(mp.dm_ind==8)); G8 = zeros(mp.Fend.corr.Npix,mp.dm8.NactTotal,mp.jac.Nmode); end % control Jacobian for DM8
if(any(mp.dm_ind==9)); G9 = zeros(mp.Fend.corr.Npix,mp.dm9.NactTotal,mp.jac.Nmode); end % control Jacobian for DM9

for tsi=1:mp.jac.Nmode
    
    Ein = ones(mp.P1.compact.Narr);
    
    modvar.sbpIndex = 1;
    modvar.wpsbpIndex = 0; %--Dummy index since not needed in compact model
    modvar.whichSource = 'star';     
    
    normFac = mp.Fend.compact.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
    flagEval = false;             % flag to use a different (usually higher) resolution at final focal plane for evaluation
    
    lambda = mp.sbp_centers(modvar.sbpIndex);
    
    Eunpoked = model_compact(mp, modvar);
    EunpokedVec = Eunpoked(mp.Fend.corr.inds);
    
    %--Define what the complex-valued FPM is if the coronagraph is some type of HLC.
    switch upper(mp.coro) 
        case{'EHLC'} %--DMs, optional apodizer, extended FPM with metal and dielectric modulation and outer stop, and LS. Uses 1-part direct MFTs to/from FPM
            mp.FPM.mask = falco_gen_EHLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'compact'); %--Complex transmission map of the FPM.
        case{'HLC','APHLC'} %--DMs, optional apodizer, FPM with optional metal and dielectric modulation, and LS. Uses Babinet's principle about FPM.
            mp.FPM.mask = falco_gen_HLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'compact'); %--Complex transmission map of the FPM.
    end

    %--DM1
    fprintf('Starting Jacobian calculation with compact model for DM1...'); tic
    if(any(mp.dm_ind==1))
        whichDM = 1;
        parfor iact=1:mp.dm1.NactTotal
            G1(:,iact,tsi) = EXAMPLE_func_validate_Jacobian_with_compact_model(iact,whichDM, Vfrac, EunpokedVec, mp,  lambda, normFac, Ein); 
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)
   
    %--DM2
    fprintf('Starting Jacobian calculation with compact model for DM2...'); tic
    if(any(mp.dm_ind==2))
        whichDM = 2;
        parfor iact=1:mp.dm1.NactTotal
            G2(:,iact,tsi) = EXAMPLE_func_validate_Jacobian_with_compact_model(iact,whichDM, Vfrac, EunpokedVec, mp,  lambda, normFac, Ein);
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)

    %--DM8
    fprintf('Starting Jacobian calculation with compact model for DM8...'); tic
    if(any(mp.dm_ind==8))
        whichDM = 8;
        Vfrac = 1;
        parfor iact=1:mp.dm8.NactTotal
            G8(:,iact,tsi) = EXAMPLE_func_validate_Jacobian_with_compact_model(iact,whichDM, Vfrac, EunpokedVec, mp,  lambda, normFac, Ein);
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)

    %--DM9
    fprintf('Starting Jacobian calculation with compact model for DM9...'); tic
    if(any(mp.dm_ind==9))
        whichDM = 9;
        Vfrac = 10;
        parfor iact=1:mp.dm9.NactTotal
            G9(:,iact,tsi) = EXAMPLE_func_validate_Jacobian_with_compact_model(iact,whichDM, Vfrac, EunpokedVec, mp,  lambda, normFac, Ein);
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)
        
end

if(any(mp.dm_ind==1));  G1_compact = G1; clear G1;  end
if(any(mp.dm_ind==2));  G2_compact = G2; clear G2;  end
if(any(mp.dm_ind==8));  G8_compact = G8; clear G8;  end
if(any(mp.dm_ind==9));  G9_compact = G9; clear G9;  end

%%
%% DM1 Comparison, Compact Model
%--Compare overall Jacobian
figure(1); imagesc(abs(G1_jac)); colorbar; set(gca,'Fontsize',20);
figure(2); imagesc(abs(G1_compact)); colorbar; set(gca,'Fontsize',20);
figure(3); imagesc(abs(G1_jac-G1_compact)); colorbar; set(gca,'Fontsize',20);
figure(4); imagesc(abs(G1_jac-G1_compact)/max(abs(G1_compact(:)))); colorbar; set(gca,'Fontsize',20); %--Normalized error

%--Compare Jacobian for just one actuator
whichAct = 685;

Etemp = model_compact(mp, modvar);
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.Fend.corr.inds) = G1_jac(:,whichAct,1);
Etemp2(mp.Fend.corr.inds) = G1_compact(:,whichAct,1);

figure(11); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(12); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(13); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(14); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error

figure(18); imagesc(angle(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(19); imagesc(angle(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);

%% DM2 Comparison, Compact Model
%--Compare overall Jacobian
figure(101); imagesc(abs(G2_jac)); colorbar; set(gca,'Fontsize',20);
figure(102); imagesc(abs(G2_compact)); colorbar; set(gca,'Fontsize',20);
figure(103); imagesc(abs(G2_jac-G2_compact)); colorbar; set(gca,'Fontsize',20);

%--Compare Jacobian for just one actuator
whichAct = 685;

Etemp = model_compact(mp, modvar);
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.Fend.corr.inds) = G2_jac(:,whichAct,1);
Etemp2(mp.Fend.corr.inds) = G2_compact(:,whichAct,1);

figure(111); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(112); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(113); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(114); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error

%%
return
%% DM8 Comparison, Compact Model

jac2D = zeros(sqrt(mp.dm8.NactTotal));
jac2D(:) = sum(abs(G8_jac).^2,1);
figure(20); imagesc(abs(jac2D)); colorbar; set(gca,'Fontsize',20);

jac2Db = zeros(sqrt(mp.dm8.NactTotal));
jac2Db(:) = sum(abs(G8_compact).^2,1);
figure(30); imagesc(abs(jac2Db)); colorbar; set(gca,'Fontsize',20);

%--Compare Jacobian for just one actuator
whichAct = 14*28+14;

Etemp = model_compact(mp, modvar); %--Just used to get the right sized image
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.Fend.corr.inds) = G8_jac(:,whichAct,1);
Etemp2(mp.Fend.corr.inds) = G8_compact(:,whichAct,1);

figure(31); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(32); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(33); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(34); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error

figure(38); imagesc(angle(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(39); imagesc(angle(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);

%% DM9 Comparison, Compact Model

%--Compare Jacobian for just one actuator
whichAct = 13*54+50;

Etemp = model_compact(mp, DM, modvar); %--Just used to get the right sized image
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.Fend.corr.inds) = G9_jac(:,whichAct,1);
Etemp2(mp.Fend.corr.inds) = G9_compact(:,whichAct,1);

figure(31); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(32); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(33); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(34); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error

figure(38); imagesc(angle(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(39); imagesc(angle(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);