function [n]=ZnSe(lam)
    % https://refractiveindex.info/?shelf=main&book=ZnSe&page=Marple
    data = csvread('ZnSe.csv',1);  % skip first row
    lams = data(:,1);              % first comumn is wavelength
    ns   = data(:,2);              % second comumn is refractive index
    n    = interp1(lams,ns,lam);   % interpolate query wavelength
end