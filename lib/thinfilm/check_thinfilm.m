

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

% From PMGI_transmission_only.xlsx
lam = 600e-9;
d0 = 4*lam;
aoi = 0;
t_Ti_base = 0;
t_Ni_vec = 0;%20e-9;%95e-9;
t_PMGI_vec = 100e-9;
pol = 2;
[tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% phase = angle(tCoef) - phase0
T = abs(tCoef)^2 % Value from Bala: 0.9431006949; Value from FALCO: 0.94313153

% From PMGIon95nmPMGI_aoi10deg_T_sPol.xlsx
lam = 400e-9;
d0 = 4*lam;
aoi = 10;
t_Ti_base = 0;
t_Ni_vec = 95e-9;
t_PMGI_vec = 0;%50e-9;
pol = 0;
[tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% phase = angle(tCoef) - phase0
T = abs(tCoef)^2 % Value from Bala: 0.00087848574  Value from FALCO: 0.000878466587


% From PMGIon95nmPMGI_aoi10deg_T_pPol.xlsx
lam = 450e-9;
d0 = 4*lam;
aoi = 10;
t_Ti_base = 0;
t_Ni_vec = 95e-9;
t_PMGI_vec = 30e-9;
pol = 1;
[tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% phase = angle(tCoef) - phase0
T = abs(tCoef)^2 % Value from Bala: 0.00118382732,   Value from FALCO: 0.00118379


% From PMGIon95nmPMGI_aoi10deg_T_pPol.xlsx
lam = 550e-9;
d0 = 4*lam;
aoi = 10;
t_Ti_base = 0;
t_Ni_vec = 95e-9;
t_PMGI_vec = 600e-9;
pol = 1;
[tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
% phase = angle(tCoef) - phase0
T = abs(tCoef)^2 % Value from Bala: 0.00121675706  Value from FALCO: 0.001216750339


% % From PMGIon95nmPMGI_aoi10deg_ph_sPol.xlsx
% lam = 500e-9;
% d0 = 4*lam;
% aoi = 10;
% t_Ti_base = 0;
% t_Ni_vec = 95e-9;
% pol = 0;
% t_PMGI_A = 10e-9;
% t_PMGI_B = 100e-9;
% 
% [tCoefA, rCoefA] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_A, d0, pol);
% [tCoefB, rCoefB] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_B, d0, pol);
% phaseA = angle(tCoefA)*180/pi;
% phaseB = angle(tCoefB)*180/pi;
% phaseB - phaseA
% 
% phaseA_Bala = -92.73355396
% phaseB_Bala = 151.6808576

