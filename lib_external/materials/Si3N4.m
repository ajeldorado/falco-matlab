function [n]=Si3N4(lam)
    % https://refractiveindex.info/?shelf=main&book=Si3N4&page=Luke
    % 310nm to 5µm
    data = csvread('Si3N4.csv',1);    % skip first row
    lams = data(:,1);                 % first comumn is wavelength
    ns   = data(:,2);                 % second comumn is refractive index
    n    = interp1(lams,ns,lam);      % interpolate query wavelength
end