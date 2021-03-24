
clear all;

Nbeam = 1000;
ID = 0.36;
OD = 0.91;
wStrut = 3.2/100.;
rocFillet = 0.02;
upsampleFactor = 100;
centering = 'pixel';
lyot = falco_gen_Roman_CGI_lyot_stop_symm_fillet(Nbeam, ID, OD, wStrut, rocFillet, upsampleFactor, centering);

figure(1); imagesc(lyot); axis xy equal tight; colorbar;