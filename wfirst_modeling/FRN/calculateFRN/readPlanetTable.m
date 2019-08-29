function planet = readPlanetTable(planetNo, bandNo)
% Reads PlanetTable.csv and returns struct of data for planet selected
%
% S. Miller 25-Jan-2019

% get path of subdirectory containing FRN scenario data
filename = mfilename('fullpath');
filepath = fileparts(filename);
scenario_Directory = fullfile(filepath, 'FRNscenarioData');

uc = unitsConstants;

% suppress warnings for variable name changes in tables
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');

planetTableMaster = readtable(fullfile(scenario_Directory, 'PlanetTable.csv'));

planetTable = planetTableMaster(planetNo, :);
planet = table2struct(planetTable);

if exist('bandNo', 'var')
    if planetNo == 1
        switch bandNo
            case 1
                % keep values in table
            case 3
                planet.A = planet.A * 4/3;
                planet.Rp_R_j = planet.Rp_R_j * 4/3;
            case 4
                planet.A = planet.A * 8/3;
                planet.Rp_R_j = planet.Rp_R_j * 8/3;
            otherwise
                error('bandNo = %.1d is not supported. Supported band numbers are 1, 3 and 4.', bandNo);
        end
    else
        error('Band numbers only supported for planet number 1 (EB Fiducial Planet)')
    end
end

% default to band number = 1 if none specified for EB Fiducial Planet
if planetNo == 1 && ~exist('bandNo', 'var')
    bandNo = 1;
end

planet.PhaseAng = planet.PhaseAng * uc.deg;
planet.FluxRatio = planet.Albedo * (planet.Rp_R_j*uc.jupiterRadius/(planet.A*uc.AU))^2;
planet.Sep = (planet.A*sin(planet.PhaseAng)*uc.AU) / (planet.PcDist*uc.parsec);
planet.bandNo = bandNo;

end
