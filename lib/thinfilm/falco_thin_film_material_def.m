% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function [tCoef] = falco_thin_film_material_def(lam, aoi, t_Ni, t_PMGI, pol)
%
% Calculates the thin-film complex transmission for the provided
% combinations of metal and dielectric thicknesses and list of wavelengths.
%
% INPUTS:
%   lam: Wavelength [m]
%   aoi:    Angle of incidense [deg]
%   t_Ni:   Nickel layer thickness [m]
%   t_PMGI: PMGI layer thickness [m]
%   pol: = 0 for TE(s) polarization, = 1 for TM(p) polarization, 2 for mean
%   of s and p polarizations
%
% OUTPUTS:
%   cMask(t_PMGI,t_ni): complex field transmission coeffient. Scalar,
%   complex value.
%
% REVISION HISTORY:
% Modified on 2019-01-28 by A.J. Riggs to:
%  -Allow for returning the mean transmission for different polarizations
%  -Add optional keyword input for OPD or non-OPD phase convention choice
%  -Add optional keyword input for substrate material choice.
%  -Cleaned up the code.
% Modified on 2018-05-01 by A.J. Riggs.
% Created on 2017-12-11 by Erkin Sidick.
% 1/25/2019: Erkin replaced Ni, Ti and Fused-Silica indices with Dwight's.
% -------------------------------------------------------------------------

function [tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol, varargin)


%% Optional Keyword Inputs

flagOPD = false; %--Default value for OPD phase sign convention is false.
substrate = 'FS'; % material name of the mask substrate [FS or N-BK7]

icav = 0;             % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'opd'}
        flagOPD  = true;       % Use the OPD phase sign convention.
      case {'substrate','sub','glass'}
        icav = icav + 1;
        substrate = varargin{icav}; % material name of the mask substrate [FS or N-BK7]
      otherwise
        error('falco_thin_film_material_def: Unknown keyword: %s\n', varargin{icav});
    end
end

%% Define Material Properties

lam_nm = lam * 1.0e9;    % m --> nm
lam_u = lam*1.0e6; % m --> microns
theta  = aoi*pi/180;     % deg --> rad

% ---------------------------------------------
%--Substrate properties
switch lower(substrate)
    case{'fs','fusedsilica','fused_silica'}   % Fused Silica
        % ----------- Fused Silica from Dwight Moody------------------
        lamm = [.4e-6,  .5e-6, .51e-6, .52e-6, .53e-6, .54e-6, .55e-6, .56e-6, .57e-6, .58e-6, .59e-6, .6e-6, .72e-6, .76e-6,   .8e-6,  .88e-6,  .90e-6 1.04e-6]*1d9;
        nx = [ 1.47012, 1.462, 1.462,  1.461,  1.461,  1.460,  1.460,  1.460,  1.459,  1.459,  1.458,  1.458, 1.45485, 1.45404, 1.45332, 1.45204, 1.45175, 1.44992];
        vsilica = [lamm(:) nx(:)];
        lam_silica = vsilica(:,1);  % nm
        n_silica   = vsilica(:,2);
        n_substrate    = interp1(lam_silica, n_silica, lam_nm, 'linear');
    
    case{'n-bk7','nbk7','bk7','bk-7'} % N-BK7
    
        B1 = 1.03961212;
        B2 = 0.231792344;
        B3 = 1.01046945;
        C1 = 0.00600069867;
        C2 = 0.0200179144;
        C3 = 103.560653;

        wvl_um = lam_u;
        n_substrate = sqrt(1 + (B1*(wvl_um).^2./((wvl_um).^2 - C1)) + (B2*(wvl_um).^2./((wvl_um).^2 - C2)) + (B3*(wvl_um).^2./((wvl_um).^2 - C3)));
    
end

% ---------------------------------------------
%--Dielectric properties
npmgi = 1.524 + 5.176e-03./lam_u.^2 + 2.105e-4./lam_u.^4;
Ndiel  = length(t_PMGI_vec);

% ---------------------------------------------
%--Metal layer properties
%--New logic: Titanium layer goes beneath Nickel only. Always include them
%together. Subtract off the thickness of the Ti layer from the intended Ni
%layer thickness.

Nmetal = length(t_Ni_vec);
t_Ti_vec = zeros(Nmetal,1);

for ii = 1:Nmetal
    if(t_Ni_vec(ii) > t_Ti_base) %--For thicker layers
        t_Ni_vec(ii) = t_Ni_vec(ii) - t_Ti_base;
    else %--For very thin layers.
        t_Ti_vec(ii) = t_Ni_vec(ii);
        t_Ni_vec(ii) = 0;
    end
end
% % GUIDE:
% if(t_Ni > t_Ti) %--For thicker layers
%     t_Ni = t_Ni - t_Ti;
% else %--For very thin layers.
%     t_Ti = t_Ni;
%     t_Ni = 0;

% from D Moody
% vnickel =...
%     [400          1.61          2.36
%     440          1.62       2.59353
%     480       1.66163       2.82958
%     500         1.678         2.966
%     510         1.697         3.023
%     520         1.716          3.08
%     530         1.735         3.137
%     540         1.754         3.194
%     550         1.773         3.251
%     560         1.792         3.308
%     570         1.811         3.365
%     580          1.83         3.423
%     590         1.849          3.48
%     600         1.869         3.537
%     640       1.98941       3.75882
%     680       2.11158       3.95737
%     720          2.25         4.115
%     760       2.38625        4.2725
%     800          2.48          4.38
%     880       2.63839       4.61452
%     900        2.6675       4.67375
%     920       2.69278       4.73667
%     960       2.74947       4.86895
%     1000       2.80976       4.99537
%     1040       2.85933       5.12178];

vnickel = load('nickel_data_from_Palik_via_Bala_wvlNM_n_k.txt');

lam_nickel = vnickel(:,1);  % nm
n_nickel   = vnickel(:,2);
k_nickel   = vnickel(:,3);
nnickel    = interp1(lam_nickel, n_nickel, lam_nm, 'linear');
knickel    = interp1(lam_nickel, k_nickel, lam_nm, 'linear');

% ---------------------------------------------
% from D Moody
titanium =[...  
    397          2.08          2.95
    413          2.14          2.98
    431          2.21          3.01
    451          2.27          3.04
    471           2.3           3.1
    496          2.36          3.19
    521          2.44           3.2
    549          2.54          3.43
    582           2.6          3.58
    617          2.67          3.74
    659          2.76          3.84
    704          2.86          3.96
    756             3          4.01
    821          3.21          4.01
    892          3.29          3.96
    984          3.35          3.97
    1088           3.5          4.02
    1216          3.62          4.15];
    
lam_ti = titanium(:,1);  % nm
n_ti   = titanium(:,2);
k_ti   = titanium(:,3);
nti    = interp1(lam_ti, n_ti, lam_nm, 'linear');
kti    = interp1(lam_ti, k_ti, lam_nm, 'linear');
% ---------------------------------------------


%% Compute the complex transmission
tCoef = zeros(Ndiel,Nmetal); %--initialize
rCoef = zeros(Ndiel,Nmetal); %--initialize
for jj = 1:Ndiel
    dpm = t_PMGI_vec(jj);
    
    for ii = 1:Nmetal
        dni = t_Ni_vec(ii);
        dti = t_Ti_vec(ii);
        
        nvec = [1 1 npmgi nnickel-1i*knickel nti-1i*kti n_substrate];
        dvec = [d0-dpm-dni-dti dpm dni dti];
        
        %--Choose polarization
        if(pol==2) %--Mean of the two
            [~, ~, rr0, tt0] = falco_thin_film_solver(nvec, dvec, theta, lam, 0);
            [~, ~, rr1, tt1] = falco_thin_film_solver(nvec, dvec, theta, lam, 1);
            rr = (rr0+rr1)/2.;
            tt = (tt0+tt1)/2.;
        elseif(pol==0 || pol==1)
            [~, ~, rr, tt] = falco_thin_film_solver(nvec, dvec, theta, lam, pol);
        else
            error('falco_thin_film_material_def.m: Wrong input value for polarization.')
        end
        
        %--Choose phase convention
        if(flagOPD==false)
            tCoef(jj,ii) = conj(tt); %--Complex field transmission coeffient, changed by erkin
            rCoef(jj,ii) = conj(rr); %--Complex field reflection coeffient, changed by erkin
        else %--OPD phase convention is negative of oppositive convention
            tCoef(jj,ii) = tt; %--Complex field transmission coeffient, changed by erkin
            rCoef(jj,ii) = rr; %--Complex field reflection coeffient, changed by erkin
        end
        
    end
end

        
end %--END OF FUNCTION

