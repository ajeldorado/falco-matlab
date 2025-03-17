function [n]=diamond(lam)
    % https://refractiveindex.info/?shelf=main&book=C&page=Phillip
    % no k beyond 220nm
    data = csvread('diamond.csv',1);  % skip first row
    lams = data(:,1);                 % first comumn is wavelength
    ns   = data(:,2);                 % second comumn is refractive index
    n    = interp1(lams,ns,lam);      % interpolate query wavelength
end