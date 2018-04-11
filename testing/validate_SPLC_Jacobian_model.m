% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to compare the SPLC Jacobian from model_Jacobian with those made by
%  differencing the outputs from model_compact and model_full.
% 
% Created by A.J. Riggs on 2018-03-23.

%% Go the correct starting directory and add all of FALCO to the Matlab path
clear;

if(~isdeployed)
  pwd0 = fileparts(which(mfilename)); %--Path to this file
  cd(pwd0);
  cd ../
  addpath(genpath(pwd)) %--To find apodizer masks and saved pupils
end


%% Step 1: Define any variable values that will overwrite the defaults (in falco_config_defaults_SPLC)

%--CRITICAL FOR JACOBIAN COMPARISON WITH FULL MODEL--Must use same F4 resolution.
%%--Final Focal Plane (F4) Properties
mp.F4.compact.res = 4; %--Pixels per lambda_c/D
mp.F4.full.res = mp.F4.compact.res; %--Pixels per lambda_c/D

% %%--Record Keeping
% mp.TrialNum = 1; %--Always use a diffrent Trial # for different calls of FALCO.
% mp.SeriesNum = 1; %--Use the same Series # for sets of similar trials.

%%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 3; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.useGPU = false; %--whether to use GPUs for Jacobian calculation

mp.controller = 'EFC';%'conEFC';  % Controller options: 'EFC' for EFC as an empirical grid search over tuning parametrs, 'conEFC' for constrained EFC using CVX.
mp.centering = 'pixel'; %--Centering on the arrays at each plane: pixel or interpixel

%%--Coronagraph and Pupil Type
mp.coro = 'SPLC';   %--Tested Options: 'Vortex','LC','SPLC'
mp.whichPupil = 'LUVOIRA5'; %--Tested options: 'WFIRST_onaxis', 'WFIRST20180103','LUVOIRA5'
mp.flagApod = true;
mp.SPname = 'luvoirA5bw10';  %--Shaped pupil identifier. Current options: '32WA194','31WA220', 'luvoirA5bw10'

%%--Pupil Plane and DM Plane Properties
mp.d_P2_dm1 = 0; % distance (along +z axis) from P2 pupil to DM1 (meters)
mp.d_dm1_dm2 = 3; % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
mp.lambda0 = 500e-9; % central wavelength of bandpass (meters)
mp.fracBW = 0.01;  % fractional bandwidth of correction (Delta lambda / lambda)
mp.Nsbp = 1; % number of sub-bandpasses across correction band 
mp.Nwpsbp = 1;% number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.


%%--Pupil Masks
switch mp.whichPupil
    case 'Simple' % Can be used to create circular and annular apertures with radial spiders 
        
        mp.P1.D = 4; %--meters, diameter of telescope (This is like HabEx A)
        mp.P1.full.Nbeam = 250; 
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 
        mp.P1.compact.Nbeam = 250;
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex. 
        
        mp.P1.IDnorm = 0;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P1.ODnorm = 1;% Outer diameter (fraction of Nbeam) 
        
        mp.P4.IDnorm = 0;% Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        mp.P4.ODnorm = 0.95;% Outer diameter (fraction of Nbeam) 
        
        mp.P1.num_strut = 0;% Number of struts 
        mp.P1.strut_angs = [];%Array of angles of the radial struts (deg)
        mp.P1.strut_width = []; % Width of the struts (fraction of pupil diam.)
        
        mp.P4.num_strut = 0;% Number of struts 
        mp.P4.strut_angs = [];%Array of angles of the radial struts (deg)
        mp.P4.strut_width = []; % Width of the struts (fraction of pupil diam.)
      
    case{'LUVOIRA5'}  % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        mp.P1.D = 15.2; %14.9760; %--meters, circumscribing diameter of telescope (used only for mas-to-lambda/D conversion)
        mp.P1.Dfac = 15.2/13.7; %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
        
        mp.P1.full.Nbeam = 500;%1000;%1000; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P1.compact.Nbeam = 250;%400;%300;%400;
        
    case 'LUVOIR_B_offaxis' % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.
        
        mp.P1.full.Nbeam = 250; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 
        
        mp.P1.compact.Nbeam = 250;
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex.

        mp.P3.IDnorm = 0;
        mp.P3.ODnorm = 0.84;
        
        mp.P4.IDnorm = 0;
        mp.P4.ODnorm = 0.82;
        
end

%%--DMs
DM.dm_ind = [1 2]; % vector of which DMs to use for control.
mp.P2.D =     46.3e-3; % beam diameter at pupil closest to the DMs  (meters)
DM.dm1.Nact = 48; % number of actuators across DM1
DM.dm2.Nact = 48; % number of actuators across DM2
DM.dm_weights = ones(9,1);   % vector of relative weighting of DMs for EFC




% %%--Controller Settings
% 
% switch mp.controller
%     case{'EFC'} % 'EFC' = empirical grid search over both overall scaling coefficient and Lagrange multiplier
%         % Take images for different Lagrange multiplier values and overall command gains and pick the value pair that gives the best contrast
%         cp.muVec = 10.^(6:-1:1);
%         cp.dmfacVec = 1;%[0.7, 1]; %[0.5, 1, 2];
%         DM.maxAbsdV = 30; %--Max +/- delta voltage step for each actuator for DMs 1,2, and/or 3
%         
%     case{'conEFC'} %--Constrained EFC, written by He Sun of Princeton University
%         DM.dm1.dVpvMax = 40;
%         DM.dm2.dVpvMax = 40;
%         cp.dmfacVec = 1;
% end
% 
% %--Voltage range restrictions
% DM.dm1.maxAbsV = 250./2.;
% DM.dm2.maxAbsV = 250./2.;
% 
% %%--Tip/Tilt Control
% mp.NlamForTT = 1; %--Number of wavelengths to control  tip/tilt at. 0,1, 2, 3, or inf (for all)
% mp.Ntt = 1; %--Number of tip/tilt offsets, including 0 (so always set >=1). 1, 4, or 5
% mp.TToffset = 1; %--tip/tilt offset (mas)
    





%% Coronagraphic Mask Properties:

% mp.flagDM1stop = false; %--logical flag whether to include the stop at DM1 or not
% mp.flagDM2stop = false; %--logical flag whether to include the stop at DM2 or not


%% Final Focal Plane (F4) Properties


% %--Specs for Correction (Corr) region and the Scoring (Score) region.
% mp.F4.corr.Rin  = 2; %--lambda0/D, inner radius of correction region
% mp.F4.score.Rin = 2; %--Needs to be <= that of Correction mask
% mp.F4.corr.Rout  = floor(DM.dm1.Nact/2*(1-mp.fracBW/2)); %--lambda0/D, outer radius of correction region
% mp.F4.score.Rout = mp.F4.corr.Rout; %--Needs to be <= that of Correction mask
% mp.F4.corr.ang  = 180; %--degrees per side
% mp.F4.score.ang = 180; %--degrees per side
% mp.F4.sides = 'both'; %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


% %%--Final Focal Plane (F4) Properties
% mp.F4.compact.res = 4; %--Pixels per lambda_c/D
% mp.F4.full.res = 4; %--Pixels per lambda_c/D
% mp.F4.FOV = 1 + mp.F4.corr.Rout; % minimum desired field of view (along both axes) in lambda0/D


%% Tip/Tilt, Spatial, and Chromatic Weighting of the Control Jacobian  #NEWFORTIPTILT
% mp.Ntt = 1; %--Number of tip/tilt offsets, including 0 (so always set >=1). 1, 4, or 5
% mp.NlamForTT = 1; %--Number of wavelengths to compute tip/tilt at. 0,1, 2, 3, or inf (for all)
% mp.TToffset = 1; %--tip/tilt offset in mas
% 
% %--Spatial pixel weighting
% mp.WspatialDef = [mp.F4.corr.Rin, mp.F4.corr.Rin+2, 1];  %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]
% 
% %--Chromatic weighting

%% Deformable Mirror (DM) Parameters

% %--DM1 parameters
% DM.dm1.Nact = 48; % number of actuators across DM1
% DM.dm1.VtoH = 1*1e-9*ones(DM.dm1.Nact); % Gains: volts to meters in surface height;
% DM.dm1.xtilt = 0;
% DM.dm1.ytilt = 0;
% DM.dm1.zrot = 0; %--clocking angle (degrees)
% DM.dm1.xc = (DM.dm1.Nact/2 - 1/2); % x-center of DM in mm, in actuator widths
% DM.dm1.yc = (DM.dm1.Nact/2 - 1/2); % x-center of DM in mm, in actuator widths
% DM.dm1.edgeBuffer = 1; % Radius (in actuator spacings) outside of pupil to compute influence functions for.
% 
% %--DM2 parameters
% DM.dm2.Nact = 48; % number of actuators across DM1
% DM.dm2.VtoH = 1*1e-9*ones(DM.dm2.Nact); % Gains: volts to meters in surface height;
% DM.dm2.xtilt = 0;
% DM.dm2.ytilt = 0;
% DM.dm2.zrot = 0; %--clocking angle (degrees)
% DM.dm2.xc = DM.dm2.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% DM.dm2.yc = DM.dm2.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% DM.dm2.edgeBuffer = 1; % Radius (in actuator spacings) outside of pupil to compute influence functions for.
% 
% %--DM Actuator characteristics
% DM.dm1.dx_inf0 = 1e-4; % meters, sampling of the influence function;
% DM.dm1.dm_spacing = 1e-3;%0.9906e-3; % meters, pitch of DM actuators
% DM.dm1.inf0 = -1*fitsread('influence_dm5v2.fits');    %  -1*fitsread('inf64.3.fits');                              
% DM.dm2.dx_inf0 = 1e-4; % meters, sampling of the influence function;
% DM.dm2.dm_spacing = 1e-3;%0.9906e-3; % meters, pitch of DM actuators
% DM.dm2.inf0 = -1*fitsread('influence_dm5v2.fits');   







%% Part 2: Call the function to define the rest of the variables and initialize the workspace
if(exist('mp','var')==false); mp.dummy = 1; end
if(exist('cp','var')==false); cp.dummy = 1; end
if(exist('ep','var')==false); ep.dummy = 1; end
if(exist('DM','var')==false); DM.dummy = 1; end
[mp,cp,ep,DM] = falco_config_defaults_SPLC(mp,cp,ep,DM); %--Load values
keyboard
%% Part 3: Save the config file    
mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(DM.dm_ind)),'DM',num2str(DM.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.F4.corr.Rin),'_OWA',num2str(mp.F4.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];
fn_config = ['data/configs/',mp.runLabel,'.mat'];
save(fn_config)
fprintf('Saved the config file: \t%s\n',fn_config)


%% Part 4: Initialize workspace 
mp.flagPlot = false;
[mp,cp,ep,DM,folders] = falco_init_ws(fn_config,mp.flagPlot);


%% Part 5: Compute the Jacobian using model_Jacobian

modvar.flagCalcJac = true; 
modvar.wpsbpIndex = mp.wi_ref;
modvar.whichSource = 'star'; 

%--Re-initialize the Jacobian arrays to full size
G1=zeros(1,1,mp.Nttlam); G2=zeros(1,1,mp.Nttlam); G3=zeros(1,1,mp.Nttlam); G4=zeros(1,1,mp.Nttlam); G5=zeros(1,1,mp.Nttlam); G6=zeros(1,1,mp.Nttlam); G7=zeros(1,1,mp.Nttlam); G8=zeros(1,1,mp.Nttlam); G9=zeros(1,1,mp.Nttlam); %--Initialize for bookkeeping in cells later 
if(any(DM.dm_ind==1)); G1 = zeros(length(mp.F4.compact.corr.inds),DM.dm1.NactTotal,mp.Nttlam); end % control Jacobian for DM1
if(any(DM.dm_ind==2)); G2 = zeros(length(mp.F4.compact.corr.inds),DM.dm2.NactTotal,mp.Nttlam); end % control Jacobian for DM2
    
%--Compute the number of total actuators for all DMs used. 

GallCell1 = {squeeze(G1(:,:,1)),squeeze(G2(:,:,1)),squeeze(G3(:,:,1)),squeeze(G4(:,:,1)),squeeze(G5(:,:,1)),squeeze(G6(:,:,1)),squeeze(G7(:,:,1)),squeeze(G8(:,:,1)),squeeze(G9(:,:,1))}; % Create the cell array. Placeholders for non-existent Jacobians to have consistent numbering
NeleAll = 0;
NeleVec = []; %--Vector of total number of used actuators for each used DM
for ii=1:numel(DM.dm_ind)
    dm_index = DM.dm_ind(ii);
    NeleAll = NeleAll + size(GallCell1{dm_index},2);
    NeleVec = [NeleVec; size(GallCell1{dm_index},2) ];
end
clear GallCell1 %--Save RAM

%--Compute the control Jacobians for each DM
jacStruct =  model_Jacobian(mp, DM);
if(any(DM.dm_ind==1)); G1_jac = jacStruct.G1; end
if(any(DM.dm_ind==2)); G2_jac = jacStruct.G2; end
clear jacStruct  %--Save RAM


%% Part 6: Compute the Jacobian using model_compact and differencing

Vfrac = 1e-4;

%--Re-initialize the Jacobian arrays to full size
if(any(DM.dm_ind==1)); G1 = zeros(length(mp.F4.compact.corr.inds),DM.dm1.NactTotal,mp.Nttlam); end % control Jacobian for DM1
if(any(DM.dm_ind==2)); G2 = zeros(length(mp.F4.compact.corr.inds),DM.dm2.NactTotal,mp.Nttlam); end % control Jacobian for DM2
    

for tsi=1:mp.Nttlam
    
    modvar.sbpIndex = mp.Wttlam_si(tsi);
    modvar.ttIndex = mp.Wttlam_ti(tsi);
    modvar.wpsbpIndex = mp.wi_ref;
    modvar.flagCalcJac = 0; 
    modvar.whichSource = 'star';     
    %lambda = mp.sbp_center_vec(modvar.sbpIndex)*mp.lamFac_vec(modvar.wpsbpIndex);
    
    %--Reset the voltage maps to zero
    DM.dm1.V = 0*DM.dm1.V;
    DM.dm2.V = 0*DM.dm2.V;
    Eunpoked = model_compact(mp, DM, modvar);
    EunpokedVec = Eunpoked(mp.F4.compact.corr.inds);
    
    
    
    %--DM1
    fprintf('Starting Jacobian calculation with compact model for DM1...'); tic
    if(any(DM.dm_ind==1));
        whichDM = 1;
        parfor iact=1:DM.dm1.NactTotal
            G1(:,iact,tsi) = func_validate_Jacobian_with_compact_model(iact,whichDM,Vfrac,EunpokedVec,mp, DM, modvar);
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)
    
    %--DM2
    fprintf('Starting Jacobian calculation with compact model for DM2...'); tic
    if(any(DM.dm_ind==2));
        whichDM = 2;
        parfor iact=1:DM.dm2.NactTotal
            G2(:,iact,tsi) = func_validate_Jacobian_with_compact_model(iact,whichDM,Vfrac,EunpokedVec,mp, DM, modvar);
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)
        
end

G1_compact = G1;
G2_compact = G2;


 %--Reset the voltage maps to zero
DM.dm1.V = 0*DM.dm1.V;
DM.dm2.V = 0*DM.dm2.V;



%% Part 7: Compare Jacobians from the Jacobian and Compact Models

% save('~/Desktop/ws_temp.mat')

%--Compare overall Jacobian
figure(1); imagesc(abs(G1_jac)); colorbar; set(gca,'Fontsize',20);
figure(2); imagesc(abs(G1_compact)); colorbar; set(gca,'Fontsize',20);
figure(3); imagesc(abs(G1_jac-G1_compact)); colorbar; set(gca,'Fontsize',20);
figure(4); imagesc(abs(G1_jac-G1_compact)/max(abs(G1_compact(:)))); colorbar; set(gca,'Fontsize',20); %--Normalized error
% figure(5); imagesc(abs(G1_jac-G1_compact)./abs(G1_compact),[0 1]); colorbar; set(gca,'Fontsize',20);

%--Compare Jacobian for just one actuator
whichAct = 1024;

Etemp = model_compact(mp, DM, modvar);
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.F4.compact.corr.inds) = G1_jac(:,whichAct,1);
Etemp2(mp.F4.compact.corr.inds) = G1_compact(:,whichAct,1);


figure(11); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(12); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(13); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(14); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error
 
% figure(15); imagesc(angle(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(real(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(imag(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);


figure(21); imagesc(angle(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(22); imagesc(angle(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);



DM.dm1.V = 0*DM.dm1.V;
DM.dm2.V = 0*DM.dm2.V;

%%
return

%% Part 8: Compute the Jacobian using model_full and differencing
%--REQUIRES THAT THE F4 PLANE RESOLUTION IS THE SAME IN THE COMPACT AND FULL MODELS!!!

Vfrac = 1e-4;

%--Re-initialize the Jacobian arrays to full size
if(any(DM.dm_ind==1)); G1 = zeros(length(mp.F4.compact.corr.inds),DM.dm1.NactTotal,mp.Nttlam); end % control Jacobian for DM1
if(any(DM.dm_ind==2)); G2 = zeros(length(mp.F4.compact.corr.inds),DM.dm2.NactTotal,mp.Nttlam); end % control Jacobian for DM2
    

for tsi=1:mp.Nttlam
    
    modvar.sbpIndex = mp.Wttlam_si(tsi);
    modvar.ttIndex = mp.Wttlam_ti(tsi);
    modvar.wpsbpIndex = mp.wi_ref;
    modvar.flagCalcJac = 0; 
    modvar.whichSource = 'star';     
    %lambda = mp.sbp_center_vec(modvar.sbpIndex)*mp.lamFac_vec(modvar.wpsbpIndex);
    
    %--Reset the voltage maps to zero
    DM.dm1.V = 0*DM.dm1.V;
    DM.dm2.V = 0*DM.dm2.V;
    Eunpoked = model_full(mp, DM, modvar);
    EunpokedVec = Eunpoked(mp.F4.compact.corr.inds);
    
    
    
    %--DM1
    fprintf('Starting Jacobian calculation with compact model for DM1...'); tic
    if(any(DM.dm_ind==1));
        whichDM = 1;
        parfor iact=1:DM.dm1.NactTotal
            G1(:,iact,tsi) = func_validate_Jacobian_with_full_model(iact,whichDM,Vfrac,EunpokedVec,mp, DM, modvar);
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)
    
    %--DM2
    fprintf('Starting Jacobian calculation with compact model for DM2...'); tic
    if(any(DM.dm_ind==2));
        whichDM = 2;
        parfor iact=1:DM.dm2.NactTotal
            G2(:,iact,tsi) = func_validate_Jacobian_with_full_model(iact,whichDM,Vfrac,EunpokedVec,mp, DM, modvar);
        end
    end
    fprintf('done. Time = %.1f sec.\n',toc)
        
end

G1_full = G1;
G2_full = G2;
clear G1 G2


 %--Reset the voltage maps to zero
DM.dm1.V = 0*DM.dm1.V;
DM.dm2.V = 0*DM.dm2.V;



%% Part 9: Compare Jacobians from the Jacobian and Compact Models

% save('~/Desktop/ws_temp.mat')

%--Compare overall Jacobian
figure(101); imagesc(abs(G1_jac)); colorbar; set(gca,'Fontsize',20);
figure(102); imagesc(abs(G1_full)); colorbar; set(gca,'Fontsize',20);
figure(103); imagesc(abs(G1_jac-G1_full)); colorbar; set(gca,'Fontsize',20);
figure(104); imagesc(abs(G1_jac-G1_full)/max(abs(G1_full(:)))); colorbar; set(gca,'Fontsize',20); %--Normalized error
% figure(5); imagesc(abs(G1_jac-G1_full)./abs(G1_full),[0 1]); colorbar; set(gca,'Fontsize',20);

%--Compare Jacobian for just one actuator
whichAct = 1024;

Etemp = model_compact(mp, DM, modvar);
Etemp1 = 0*Etemp;
Etemp2 = 0*Etemp;
Etemp1(mp.F4.compact.corr.inds) = G1_jac(:,whichAct,1);
Etemp2(mp.F4.compact.corr.inds) = G1_full(:,whichAct,1);


figure(111); imagesc(abs(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(112); imagesc(abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(113); imagesc(abs(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(114); imagesc(abs(Etemp1-Etemp2)./abs(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20); %--Normalized error
 
% figure(15); imagesc(angle(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(real(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
% figure(16); imagesc(imag(Etemp1-Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);


figure(121); imagesc(angle(Etemp1)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);
figure(122); imagesc(angle(Etemp2)); axis xy equal tight; colorbar; set(gca,'Fontsize',20);



DM.dm1.V = 0*DM.dm1.V;
DM.dm2.V = 0*DM.dm2.V;





