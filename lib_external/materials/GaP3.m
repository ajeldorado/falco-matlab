function [n]=GaP3(lam)
    % https://refractiveindex.info/?shelf=main&book=GaP&page=Khmelevskaia
    % Gallium Phosphide. transparent beyond lam=690nm. Data up to 1700nm.
    data = csvread('GaP3.csv',1);  % skip first row
    lams = data(:,1);                 % first comumn is wavelength
    ns   = data(:,2);                 % second comumn is refractive index
    n    = interp1(lams,ns,lam);      % interpolate query wavelength
end
