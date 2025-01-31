% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform a DMVC simple design run.
%  1) Load the default model parameters for a vortex.
%  2) Specify the values to overwrite.
%  3) Run a single trial of WFC using FALCO.


clear
% close all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
%mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
%mp.path.proper = '~/Documents/MATLAB/proper_v3.2_matlab_11feb20/'; %--Location of the MATLAB PROPER library

%mp.path.falco = '~/Documents/Documents/FALCO_MATLAB/falco-matlab-master_2022';  %--Location of FALCO
%mp.path.proper = '/Applications/PROPER/'; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
%mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
%mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%mp.path.config = '~/Documents/Documents/FALCO_MATLAB/falco-matlab-master_2022/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
%mp.path.ws = '~/Documents/Documents/FALCO_MATLAB/falco-matlab-master_2022/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
%addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
%addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


%% Step 2: Load default model parameters

EXAMPLE_defaults_OMC

mp.minimizeNI= true;
mp.flagSVD = true;

%% Step 3: Overwrite default values as desired

% %%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = '/Users/cprada/Documents/Documents/FALCO_MATLAB/falco-matlab-master_2022/data/brief/Series0001_Trial0001_LC_simple_1DM48_z1_IWA2.8_OWA10_1lams550nm_BW1_plannedEFC_snippet.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.dm1.Vall(:,:,31);
% mp.dm2.V = temp.out.dm2.Vall(:,:,31);
% clear temp

% mp.jac.star.weights = [1, 1]; % star weights for control
% 


% %--Correction and scoring region definition
% mp.Fend.corr.Rin = 0;   % inner radius of dark hole correction region [lambda0/D]
% mp.Fend.corr.Rout  = 2;  % outer radius of dark hole correction region [lambda0/D]
% mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]
% mp.Fend.score.Rin = 0;  % inner radius of dark hole scoring region [lambda0/D]
% mp.Fend.score.Rout = 2;  % outer radius of dark hole scoring region [lambda0/D]
% mp.Fend.score.ang = 180;  % angular opening of dark hole scoring region [degrees]
% mp.Fend.shape = 'square';
% mp.Fend.xiOffset = 6;

mp.ctrl.log10regVec     = -6:1/2:0; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.log10regSchedIn = [-2.000  -2.000   -2.000   -2.000   -2.000   ...
                           -2.000  -2.000   -2.000   -2.000   -2.000   ...
                            1j     -2.000   -2.000   -2.000   -2.000   ...
                            1j     -2.500   -2.500   -2.500   -2.500   ...
                            1j   1j  1j   1j   1j        ...
                            1j   1j  1j   1j   1j ];
%mp.ctrl.flagUseModel = false;%true;

%%--Special Computational Settings
mp.flagParfor = true;%false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 3;

%--Use just 1 wavelength for initial debugging of code
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass


%mp.Nitr = 29; %--Number of wavefront control iterations

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control


[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);


%%%%%%%% plot the contrast

figure(991);
semilogy(0:out.Itr,out.InormHist);
xlabel('Iteration')
ylabel('Normalized Intensity')

% %%%%%%%%%   This would be how to plot the wavelengths individually
% 
%  figure(991);
% for i = 1:5
% semilogy(1:out.Itr,out.normIntMeasScore(:,i));hold on; 
% end

% %%%%% Likewise, you can look at the normalized intensity vs wavelength:
% 
% figure(992);
% semilogy(mp.sbp_centers*1e9,out.normIntMeasScore(end,:));
% xlabel('Wavelength (nm)');
% ylabel('Normalized intensity')
% 
% %%%%%%%%%%%%This is how you would plot the the regularization history
% 
figure(993);
plot(1:out.Itr,out.log10regHist)
xlabel('Iteration')
ylabel('Regularization')
