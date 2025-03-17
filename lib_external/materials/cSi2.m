function [n]=cSi(lam)
    % https://refractiveindex.info/?shelf=main&book=Si&page=Li-293K
    % alternative: Cristalline Si. n and k data starting at 1.2mu upwards
    data = csvread('cSi2.csv',1);  % skip first row
    lams = data(:,1);                 % first comumn is wavelength
    ns   = data(:,2);                 % second comumn is refractive index
    n    = interp1(lams,ns,lam);      % interpolate query wavelength
end
