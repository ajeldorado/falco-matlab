

clear all

% 
% lam = 400e-9;
% d0 = 4*lam;
% aoi = 10;
% t_Ti_base = 0;
% t_Ni_vec = 95e-9;
% t_PMGI_vec = 0;
% pol = 2;
% 
% 
% [tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol)


% lam = 400e-9;
% d0 = 4*lam;
% aoi = 10;
% t_Ti_base = 0;
% t_Ni_vec = 20e-9;%95e-9;
% t_PMGI_vec = 0;
% pol = 2;
% 
% 
% [tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol)
% 
% abs(tCoef)^2


% lam = 400e-9;
% d0 = 4*lam;
% aoi = 10;
% t_Ti_base = 0;
% t_Ni_vec = 0;%20e-9;%95e-9;
% t_PMGI_vec = 0;%40e-9;
% pol = 2;
% [tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% % phase0 = angle(tCoef)
% T = abs(tCoef)^2 % Value from Bala: 96.3687853267053

% lam = 500e-9;
% d0 = 4*lam;
% aoi = 10;
% t_Ti_base = 0;
% t_Ni_vec = 0;%20e-9;%95e-9;
% t_PMGI_vec = 50e-9;
% pol = 2;
% [tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% % phase = angle(tCoef) - phase0
% T = abs(tCoef)^2 % Value from Bala: 94.8794827030582


% lam = 600e-9;
% d0 = 4*lam;
% aoi = 10;
% t_Ti_base = 0;
% t_Ni_vec = 0;%20e-9;%95e-9;
% t_PMGI_vec = 100e-9;
% pol = 2;
% [tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% % phase = angle(tCoef) - phase0
% T = abs(tCoef)^2 % Value from Bala: 94.31006949


lam = 400e-9;
d0 = 4*lam;
aoi = 10;
t_Ti_base = 0;
t_Ni_vec = 95e-9;
t_PMGI_vec = 0;%50e-9;
pol = 2;
[tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% phase = angle(tCoef) - phase0
T = abs(tCoef)^2 % Value from Bala: 94.8794827030582



lam = 450e-9;
d0 = 4*lam;
aoi = 10;
t_Ti_base = 0;
t_Ni_vec = 95e-9;
t_PMGI_vec = 30e-9;
pol = 2;
[tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% phase = angle(tCoef) - phase0
T = abs(tCoef)^2 % Value from Bala: 0.118382731718367,   Value from FALCO: 0.1177

