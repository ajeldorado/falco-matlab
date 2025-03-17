function [n]=GaP1(lam)
    % https://refractiveindex.info/?shelf=main&book=GaP&page=Aspnes
    % Gallium Phosphide. transparent beyond lam=564nm. Data up to 826nm.
    data = csvread('GaP1.csv',1);  % skip first row
    lams = data(:,1);                 % first comumn is wavelength
    ns   = data(:,2);                 % second comumn is refractive index
    n    = interp1(lams,ns,lam);      % interpolate query wavelength
end
