function [n]=GaP2(lam)
    % https://refractiveindex.info/?shelf=main&book=GaP&page=Jellison
    % Gallium Phosphide. transparent in range lam=480-825nm. Data up to 840nm.
    data = csvread('GaP2.csv',1);  % skip first row
    lams = data(:,1);                 % first comumn is wavelength
    ns   = data(:,2);                 % second comumn is refractive index
    n    = interp1(lams,ns,lam);      % interpolate query wavelength
end
