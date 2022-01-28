%% Define Necessary Paths on Your System

%--Library locations
mp.path.falco = '/Users/jllopsay/Documents/GitHub/falco-matlab/';  %--Location of FALCO
mp.path.proper = '/Users/jllopsay/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library

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
fn_config = '/Users/jllopsay/Documents/GitHub/falco-matlab/data/brief/Series0002_Trial0001_vortex_Simple_1DM34_z0.2_IWA3_OWA10_1lams775nm_BW1_gridsearchEFC_config.mat';
fn_snippet = '/Users/jllopsay/Documents/GitHub/falco-matlab/data/brief/Series0002_Trial0001_vortex_Simple_1DM34_z0.2_IWA3_OWA10_1lams775nm_BW1_gridsearchEFC_snippet.mat';

%% Compute sensitivities to 1nm RMS of the specified Zernikes
dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet);
% dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet,'plot');
% dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet,'parfor');
% dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet,'plot','parfor');