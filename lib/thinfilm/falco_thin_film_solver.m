% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute thin film equations with 2x2 transfer matrices.
%
% [R, T, rr, tt] = falco_thin_film_solver(n,d,theta,lam,[tetm])

% INPUTS
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
% OUTPUTS
% R = normalized reflected intensity coefficient
% T = normalized transmitted intensity coefficient
% rr = complex field reflection coefficient
% tt = complex field transmission coefficient

function [R, T, rr, tt] = falco_thin_film_solver(n,d,theta,lam,tetm)

    % Validate inputs
    if nargin == 0
        error('usage: [R, T, rr, tt] = falco_thin_film_solver(n,d,theta,lam,[tetm])');
    end

    N = length(n);
    if length(d) ~= N-2
        error('n and d mismatch');
    end
    n = n(:); d = [0; d(:); 0];

    if nargin < 5, tetm = 0; end

    % calculate wave vector of incident wave, 
    %    kx.^2 + kz.^2 = k.^2 = (2 pi n_incident_medium / wavelength)^2
    kx = 2*pi*n(1)*sin(theta)/lam;
    kz = -sqrt( (2*pi*n/lam).^2 - kx.^2 ); % sign agrees with measurement convention

    % for TE, E = Ey (perpendicular to the plane of incidence). Boundary condition
    %         at each layer is continuity of tangential E = Ey(z=constant), and 
    %         tangential H = d/dz (Ey) = -kz Ey(z=constant) (assumes all materials are
    %         non-magnetic, mu = mu_0 vacuum permeability)
    %         Since x is not constant along the boundary, kx is constant through the
    %         the whole stack, and kz is calculated for each layer (=kzz below)
    % for TM, H = Hy (perpendicular to the plane of incidence). Boundary condition
    %         at each layer is continuity of tangential H = Hy, and
    %         tangential E. From Maxwell's eqn:
    %         eps E = curl H
    %         E_tangential = (1/eps) d/dz Hy = (1/n^2)kz Hy
    %         Therefore, kzz = kz./(n.^2) for each layer (n is for each
    %         layer)
    if tetm == 1
       kzz = kz./(n.^2);
    else
       kzz = kz;
    end

    % in each layer, there is a forward propagating wave, Et, and a backward
    % propagating wave Er:
    eep = exp(-1i*kz.*d);
    eem = exp(1i*kz.*d);

    % Snell's law from the boundary conditions give a r and t coefficient at
    % each boundary. Note that kzz is appropriate for either TE or TM fields
    % from above.
    i1  = 1:N-1; 
    i2  = 2:N;
    tin = 0.5*(kzz(i1) + kzz(i2))./kzz(i1); % Snell's law for transmission
    ri  = (kzz(i1) - kzz(i2))./(kzz(i1) + kzz(i2)); % Snell's law for refl

    % Here's the tricky part. At each boundary, the boundary condition is:
    % Ft(i) + Fr(i) = Ft(i+1) + Fr(i+1),   (eqn 1)
    % where F is either Ey for TE or Hy for TM
    % t = transmitted forward propagating wave,
    % r = reflected backward propagating wave
    % F(i)   is the field in the medium to the left of the boundary
    % F(i+1) is the field in the medium to the right of the boundary
    % Then, coupled equations at the boundary:
    % eem * Fr(i) = eep * ri *Ft(i) + tin *Fr(i+1)
    % Ft(i+1) = eep * tin * Ft(i) - ri * Fr(i+1)
    % where we chose a convention that the fields on the left size have zero
    % phase, and the fields on the right side are phase retarded by the layer
    % to the right. Note that the transmission coefficient is the same for
    % left-to-right or right-to-left, but the reflection coefficient is
    % opposite sign for right-to-left than for left-to-right.
    % Solve the coupled equations to make a matrix equation:
    % [Et(i); Er(i)] = A * [Et(i+1); Er(i+1)]
    % matrix multiply the A for each boundary:
    A = eye(2);
    for i = 1:N-1
        A = A * (tin(i)*[eep(i) ri(i)*eep(i); ri(i)*eem(i) eem(i)]);
    end

    % now have [Ft(0); Fr(0)] = A * [Ft(N); Fr(N)]
    % but Ft(0) = incident field = 1,
    % and Fr(N) = 0.
    % then
    rr = A(2, 1)/A(1, 1); % = Fr(0) = reflection coefficent in incident
    tt = 1/A(1, 1); % = Ft(N) = transmission coefficient in substrate

    % transmitted power flux (Poynting vector . surface) depends on index of the
    % substrate and angle
    R = abs(rr).^2;
    if tetm == 1
        Pn = real( (kz(N)/(n(N).^2)) ./ (kz(1)/(n(1).^2)) );
    else
        Pn = real((kz(N)./kz(1)));
    end
    T = Pn.*abs(tt).^2;
    tt= sqrt(Pn).*tt;

return
