% -------------------------------------------------------------------------
%
% Script to return index of refraction for a material at a certain wavelength


% REVISION HISTORY:
% --------------
% Created on 2025-1-3 by Niyati Desai.

% lam input should be in microns
% Material options are (included in /lib_external/materials folder):
% TiO2
% aSi_k
% aSi (alpha silicon)
% cSi (crystalline silicon)
% diamond
% GaP1
% HfO2
% Si3N4
% SiO2
% TiO2
% ZnSe
% fusedsilica


function [n]=getn(lam,material)    
    csvfile = material+".csv";
    
    data = csvread(csvfile,1);  % skip first row
    lams = data(:,1);              % first comumn is wavelength
    ns   = data(:,2);              % second comumn is refractive index
    n    = interp1(lams,ns,lam);   % interpolate query wavelength
end