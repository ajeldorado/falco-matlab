function [n]=SiO2(lam)
    % https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson
    % fused silica (SiO2). no k.
    data = csvread('SiO2.csv',1);  % skip first row
    lams = data(:,1);                 % first comumn is wavelength
    ns   = data(:,2);                 % second comumn is refractive index
    n    = interp1(lams,ns,lam);      % interpolate query wavelength
end
