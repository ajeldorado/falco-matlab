clear

Nbeam = 500;
ID = 0.50;
OD = 0.80;
wStrut = 3.6/100.;
rocFillet = 0.02;
upsampleFactor = 100;
centering = 'pixel';
lyotRounded = falco_gen_Roman_CGI_lyot_stop_symm_fillet(Nbeam, ID, OD, wStrut, rocFillet, upsampleFactor, centering);

changes.flagLyot = true;
changes.ID = ID;
changes.OD = OD;
changes.wStrut = wStrut;
lyotSharp = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);

nArray = ceil_odd(length(lyotRounded));
lyotRounded = pad_crop(lyotRounded, nArray);
lyotSharp = pad_crop(lyotSharp, nArray);

figure(1); imagesc(lyotRounded); axis xy equal tight; colorbar;
figure(2); imagesc(lyotSharp); axis xy equal tight; colorbar;
figure(3); imagesc(lyotSharp - lyotRounded); axis xy equal tight; colorbar;
