function uc = unitsConstants 

%-------------------------------------------------------------------------
% Units 
% All quantities are converted and held in SI units, and angles in radians

uc.ppb = 1e-9; % parts per billion

uc.meter = 1; uc.second = 1; uc.kg = 1;
uc.C = 1; uc.Ohm = 1; uc.Farad = 1; uc.Kelvin = 1; uc.Joule = 1; uc.Watt = 1; 
uc.radian = 1;
uc.km = 1000 * uc.meter; uc.cm = .01 * uc.meter; uc.mm = .001* uc.meter; uc.inch = 2.54*uc.cm;
uc.um = 1e-6*uc.meter; uc.nm = 1e-9*uc.meter;
uc.mW = 1e-3*uc.Watt; uc.uW = 1e-6*uc.Watt; uc.nW = 1e-9*uc.Watt;
uc.Ampere = 1 * uc.C/uc.second; uc.mA = 1e-3 * uc.Ampere;  uc.uA = 1e-6 * uc.Ampere; 

uc.Hz = 1/uc.second; uc.kHz = 1e3 * uc.Hz;
uc.deg = (pi / 180) * uc.radian ; uc.mrad = 0.001 * uc.radian; uc.urad = 0.001 * uc.mrad;
uc.arcsec = uc.deg / 3600; uc.mas = uc.arcsec / 1000; uc.uas = uc.mas / 1000;

uc.hour = 3600 * uc.second;
uc.day  = 24 * uc.hour;
uc.usec = 1e-6*uc.second;

uc.h_planck    = 6.62607e-34; % * meter^2 * kg / second;
uc.c_light     = 299792458; % * meter / second;
uc.jupiterRadius = 69911000; % * meter;
uc.earthRadius   = 6371000; % * meter;
uc.sunAbsMag = 4.83;
uc.lightyear = 9.4607e15; % * meter;
uc.AU = 149597870700 * uc.meter;
uc.parsec = 30855152366503100;

return
