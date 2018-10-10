% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% thin_film_filter_2.m -> falco_thin_film_filter.m
%
% Function from David Marx to compute thin film equations.


function [R, T, rr, tt] = falco_thin_film_solver(n,d,theta,lam,tetm)
% [R, T, rr, tt] = thin_film_filter_2(n,d,theta,lam,[tetm])
% n = index of refraction for each layer. 
%     n(1) = index of incident medium
%     n(N) = index of transmission medium
%     then length(n) must be >= 2
% d = thickness of each layer, not counting incident medium or transmission
%     medium. length(d) = length(n)-2
% theta = angle of incidence [rad], scalar only
% lam = wavelength (scalar only). units of lam must be same as d
% tetm: 0 => TE (default), 1 => TM
%
% outputs:
% R = normalized reflected intensity coefficient
% T =        "   transmitted     "
% rr = complex field reflection coefficient
% tt =      "        transmission "

if nargin == 0, error('usage: [R, T, rr, tt] = falco_thin_film_solver(n,d,theta,lam,[tetm])'); end

N = length(n);
if length(d) ~= N-2, error('n and d mismatch'); end
n = n(:); d = [0; d(:); 0];

if nargin < 5, tetm = 0; end

kx = 2*pi*n(1)*sin(theta)/lam;
kz = -sqrt( (2*pi*n/lam).^2 - kx.^2 ); % sign agrees with measurement convention

if tetm == 1,
   kzz = kz./(n.^2);
else,
   kzz = kz;
end

eep = exp(-1i*kz.*d);
eem = exp(1i*kz.*d);

i1  = 1:N-1; 
i2  = 2:N;
tin = 0.5*(kzz(i1) + kzz(i2))./kzz(i1);
ri  = (kzz(i1) - kzz(i2))./(kzz(i1) + kzz(i2));

A = eye(2);
for i = 1:N-1,
	A = A * (tin(i)*[eep(i) ri(i)*eep(i); ri(i)*eem(i) eem(i)]);
end

rr = A(2,1)/A(1,1);
tt = 1/A(1,1);

% transmitted power flux (Poynting vector . surface) depends on index of the
% substrate and angle
R = abs(rr).^2;
if tetm == 1,
	Pn = real( (kz(N)/(n(N).^2)) ./ (kz(1)/(n(1).^2)) );
else,
	Pn = real((kz(N)./kz(1)));
end
T = Pn.*abs(tt).^2;
tt= sqrt(Pn).*tt;

return

%%%%%%%%% Old method
A = eye(2);
for i = 1:N-1,
   a = [ eem(i) eep(i); -kzz(i)*eem(i) kzz(i)*eep(i) ];
   A = A * inv(a) * ([ 1 1; -kzz(i+1) kzz(i+1)]);
end
