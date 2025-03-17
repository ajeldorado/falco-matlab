function [n]=cSi(lam)
    % https://refractiveindex.info/?shelf=main&book=Si&page=Schinke
    % Cristalline Si. n and k data up to only lam=1.45
    data = csvread('cSi.csv',1);  % skip first row
    lams = data(:,1);                 % first comumn is wavelength
    ns   = data(:,2);                 % second comumn is refractive index
    n    = interp1(lams,ns,lam);      % interpolate query wavelength
end