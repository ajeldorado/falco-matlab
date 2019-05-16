% Function to generate the polarization aberrations for Cycle 6 of the
% WFIRST CGI in the Phase A model.
%
% Translated by A.J. Riggs on 2017-06-26 from the IDL code polab.pro by John Krist.

% pro polab, polfile, lambda_m, pupil_diam_pix, condition, amp, pha
% 
% ;-- polfile: rootname of file containing polarization coefficients
% ;-- lambda_m: wavelength in meters
% ;-- pupil_diam_pix: diameter of pupil in pixels
% ;-- condition: polarization circumstance:
% ;--		-2: -45 deg in, Y out
% ;--		-1: -45 deg in, X out
% ;--		 1: +45 deg in, X out
% ;--		 2: +45 deg in, Y out
% ;-- amp, pha: returned aberration maps (pha is WFE in meters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [amp,pha] = wfirst_polab(polfile, lambda_m, pupil_diam_pix, condition)

% %--Temporarily hard-code input values in order to test as a script first.
% polfile = 'cycle6_polarization_'; % rootname of file containing polarization coefficients
% lambda_m = 910e-9;
% pupil_diam_pix = 250;
% condition = -2;

%--Choose 
if ( abs(condition) == 1 );  dir_out = 1; else dir_out = 2; end
if ( condition < 0 );  dir_in = 1; else dir_in = 2; end

%--Read in known points for polarization aberrations
zamp_array = fitsread([polfile 'amp.fits']);
zpha_array = fitsread([polfile 'pha.fits']);
lam_array_m = (450:100:950)*1e-9;

%--Interpolate to get zernike coefficients at the specified wavelength
zamp = zeros(22,1);
zpha = zeros(22,1);
for iz = 1:22 
    zamp(iz) = spline( lam_array_m, squeeze(zamp_array(:,iz,dir_in,dir_out)), lambda_m );
    zpha(iz) = spline( lam_array_m, squeeze(zpha_array(:,iz,dir_in,dir_out)), lambda_m );
end

%--Define output array size and coordinates
n = 2*floor(1/2 * pupil_diam_pix * 1.1);
x_vec = ( (0:(n-1)) - n/2 ) / (pupil_diam_pix / 2.);

%--Compute the amplitude and phase maps as matrices directly
amp = zeros(n,n);
pha = zeros(n,n);

[x,y] = meshgrid(x_vec);
r2 = x.^2 + y.^2;
r = sqrt(r2);
r3 = r.^3;
r4 = r.^4;
r5 = r.^5;
r6 = r.^6;
t = atan2(y,x);

%--Sum up the zernikes
for itype = 0:1 		%-- 0 = amp, 1 = phase
    map = zeros(n);%zeros(1,n);

    if( itype == 0 )
        z = zamp;
        map = map + z(1);	%-- include piston if amplitude map
    else
        z = zpha;
    end

    map = map + z(2) * 2 * x;				%-- x tilt
    map = map + z(3) * 2 * y;				%-- y tilt
    map = map + z(4) * sqrt(3) * (2*r2 - 1);			%-- focus
    map = map + z(5) * sqrt(6) .* r2 .* sin(2*t);		%-- 45 deg astig
    map = map + z(6) * sqrt(6) .* r2 .* cos(2*t);		%-- 0 deg astig
    map = map + z(7) * sqrt(8) .* (3*r3 - 2*r) .* sin(t);  	%-- y coma
    map = map + z(8) * sqrt(8) .* (3*r3 - 2*r) .* cos(t);	%-- x coma
    map = map + z(9) * sqrt(8) .* r3 .* sin(3*t);		%-- y trefoil 
    map = map + z(10) * sqrt(8) .* r3 .* cos(3*t);		%-- x trefoil 
    map = map + z(11) * sqrt(5) .* (6*r4 - 6*r2 + 1);		%-- spherical
        map = map + z(12) * sqrt(10) .* (4*r4 - 3*r2) .* cos(2*t);
        map = map + z(13) * sqrt(10) .* (4*r4 - 3*r2) .* sin(2*t);
        map = map + z(14) * sqrt(10) .* r4 .* cos(4*t);
        map = map + z(15) * sqrt(10) .* r4 .* sin(4*t);
        map = map + z(16) * sqrt(12) .* (10*r5 - 12*r3 + 3*r) .* cos(t);
        map = map + z(17) * sqrt(12) .* (10*r5 - 12*r3 + 3*r) .* sin(t);
        map = map + z(18) * sqrt(12) .* (5*r5 - 4*r3) .* cos(3*t);
        map = map + z(19) * sqrt(12) .* (5*r5 - 4*r3) .* sin(3*t);
        map = map + z(20) * sqrt(12) .* r5 .* cos(5*t);
        map = map + z(21) * sqrt(12) .* r5 .* sin(5*t);
        map = map + z(22) * sqrt(7) .* (20*r6 - 30*r4 + 12*r2 - 1);

    if( itype == 0 ); 
        amp = map; 
    else
        pha = map; 
    end

end

%--Plots for debugging
% mask = ones(size(r));
% mask(r>1) = 0;
% figure(101); imagesc(mask.*amp); axis xy equal tight; colorbar;
% figure(102); imagesc(x_vec,x_vec,mask.*pha/lambda_m); axis xy equal tight; colorbar;

% %--Compute the maps row by row (the way it was done in IDL)
% for jj = 1:n 
% 	y = x(jj);
% 	r2 = x.^2 + y.^2;
% 	r = sqrt(r2);
% 	r3 = r.^3;
% 	r4 = r.^4;
% 	r5 = r.^5;
% 	r6 = r.^6;
% 	t = atan2(y,x);
% 
%     for itype = 0:1 		%-- 0 = amp, 1 = phase
% 		map = zeros(1,n);
% 
%         if( itype == 0 )
% 			z = zamp;
% 			map = map + z(1);	%-- include piston if amplitude map
%         else
% 			z = zpha;
%         end
%         
% 		map = map + z(2) * 2 * x;				%-- x tilt
% 		map = map + z(3) * 2 * y;				%-- y tilt
% 		map = map + z(4) * sqrt(3) * (2*r2 - 1);			%-- focus
% 		map = map + z(5) * sqrt(6) .* r2 .* sin(2*t);		%-- 45 deg astig
% 		map = map + z(6) * sqrt(6) .* r2 .* cos(2*t);		%-- 0 deg astig
% 		map = map + z(7) * sqrt(8) .* (3*r3 - 2*r) .* sin(t);  	%-- y coma
% 		map = map + z(8) * sqrt(8) .* (3*r3 - 2*r) .* cos(t);	%-- x coma
% 		map = map + z(9) * sqrt(8) .* r3 .* sin(3*t);		%-- y trefoil 
% 		map = map + z(10) * sqrt(8) .* r3 .* cos(3*t);		%-- x trefoil 
% 		map = map + z(11) * sqrt(5) .* (6*r4 - 6*r2 + 1);		%-- spherical
%       		map = map + z(12) * sqrt(10) .* (4*r4 - 3*r2) .* cos(2*t);
%       		map = map + z(13) * sqrt(10) .* (4*r4 - 3*r2) .* sin(2*t);
%       		map = map + z(14) * sqrt(10) .* r4 .* cos(4*t);
%       		map = map + z(15) * sqrt(10) .* r4 .* sin(4*t);
%       		map = map + z(16) * sqrt(12) .* (10*r5 - 12*r3 + 3*r) .* cos(t);
%       		map = map + z(17) * sqrt(12) .* (10*r5 - 12*r3 + 3*r) .* sin(t);
%       		map = map + z(18) * sqrt(12) .* (5*r5 - 4*r3) .* cos(3*t);
%       		map = map + z(19) * sqrt(12) .* (5*r5 - 4*r3) .* sin(3*t);
%       		map = map + z(20) * sqrt(12) .* r5 .* cos(5*t);
%       		map = map + z(21) * sqrt(12) .* r5 .* sin(5*t);
%       		map = map + z(22) * sqrt(7) .* (20*r6 - 30*r4 + 12*r2 - 1);
% 
%         if( itype == 0 ); 
%             amp(jj,:) = map; 
%         else
%             pha(jj,:) = map; 
%         end
%         
%     end
% end



