%% Define Necessary Paths on Your System

%--Library locations
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


%%   ***   THE ONLY VARIABLES YOU NEED TO CHANGE ARE IN THIS CELL    ***

%--Specify zernike modes and their RMS values to use.
indsZnoll = 2:6; %--Noll indices of Zernikes to compute values for

%--Annuli to compute Zernike sensitivities over. 
% Columns are [inner radius, outer radius]. One row per annulus.
Rsens = [3, 4;...
         4, 8];

%--Input files from the FALCO trial (the configuration data and final result snippet)
fn_config = '~/Repos/falco-matlab/data/brief/Series0029_Trial0004_HLC_WFIRST180718_3DM48_z1_IWA2.7_OWA10_5lams575nm_BW10_plannedEFC_config.mat';
fn_snippet = '~/Repos/falco-matlab/data/brief/Series0029_Trial0004_HLC_WFIRST180718_3DM48_z1_IWA2.7_OWA10_5lams575nm_BW10_plannedEFC_snippet.mat';


%% Compute sensitivities to 1nm RMS of the specified Zernikes
dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet);
% dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet,'plot');
% dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet,'parfor');
% dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet,'plot','parfor');