function [table, data, QEcurves, strayLightT] = readDataTables(CG_Directory, MUFcaseInput, mode, centerLambda, bandWidth, planetWA, uc)
% Read scenario and coronagraph files into tables and extract data
% 
% S. Miller 28-Jan-2019

% get path of subdirectory containing FRN scenario data
filename = mfilename('fullpath');
filepath = fileparts(filename);
scenario_Directory = fullfile(filepath, 'FRNscenarioData');

% suppress warnings for variable name changes in tables
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');

switch MUFcaseInput
    case 'Standard'
        MUFcase = 0;
    case 'Unity'
        MUFcase = 1;
    otherwise
        error('MUFcaseInput: Please choose either Unity or Standard')
end

%% Read from FRN_scenario_data directory
% Disturbance Table
distT = readtable(fullfile(scenario_Directory, 'Disturbance.csv'), 'ReadRowNames',true);
distUnitsT = readtable(fullfile(scenario_Directory, 'DisturbanceUnits.csv'));
dist_tmp = table2array(distT);
dist = dist_tmp(dist_tmp(:,4)==MUFcase,5);
distUnits = table2array(distUnitsT);

% HST Synphot Spectra
[HSTsynT, spMag0] = readHSTsynPhotSpectra(scenario_Directory, 'HSTsynPhotSpectra.csv', centerLambda, bandWidth, uc);

% Throughput Table
[ThroughputT, t_refl] = readThroughputTable(scenario_Directory, 'ThroughputTable.csv', mode, centerLambda);

% Stray Light Table
strayLightT = readtable(fullfile(scenario_Directory, 'StrayLight.csv'), 'ReadRowNames',true);

%% Read from CG_data directory
% Lookup annular zone based on planetWA
annzoneListT = readtable(fullfile(CG_Directory, 'AnnZoneList.csv'));
annzoneList  = table2array(annzoneListT);

planetWAbounded = min(annzoneList(5,3), max(annzoneList(1,3),floor(planetWA)));
planetWAListMin = annzoneList(1:5, 3);
annzone = annzoneList(planetWAbounded == planetWAListMin);

% Initial Raw Contrast Table
initRCT = readtable(fullfile(CG_Directory, 'InitialRawContrast.csv'), 'ReadRowNames',true);
initRC_tmp = table2array(initRCT);
initRC = initRC_tmp((initRC_tmp(:,2)==annzone) & (initRC_tmp(:,3)==MUFcase), 4);

% Krist Table
kristT = readDesignFile(CG_Directory, 'KristTable.csv');

% Sensitivity Table
sensT = readtable(fullfile(CG_Directory, 'Sensitivities.csv'), 'ReadRowNames',true);
sensMUFT = readtable(fullfile(CG_Directory, 'SensitivityMUF.csv'), 'ReadRowNames',true);
sens_tmp = table2array(sensT);
sens = sens_tmp(sens_tmp(:,2)==annzone,3) .* table2array(sensMUFT(:,MUFcase+2));

%% Read QE table
QEcurvesTable = load(fullfile(scenario_Directory, 'QE_e2v.mat'), 'QEcurves');

%% Output
table.scenario.Disturbance = distT;
table.scenario.DistubanceUnits = distUnitsT;
table.scenario.HSTsynPhotSpectra = HSTsynT;
table.scenario.ThroughputTable = ThroughputT;
table.CG.AnnZoneList = annzoneListT;
table.CG.InitialRawContrast = initRCT;
table.CG.KristTable = kristT;
table.CG.Sensitivities = sensT;
table.CG.SensitivityMUF = sensMUFT;

data.annzone = annzone;
data.t_refl = t_refl;
data.spMag0 = spMag0;
data.dist = dist;
data.distUnits = distUnits;
data.initRC = initRC;
data.sens = sens;

QEcurves = QEcurvesTable.QEcurves;

return
